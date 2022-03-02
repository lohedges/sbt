/*
  Copyright (c) 2021-2022 Lester Hedges <lester.hedges@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "EventReader.h"
#include "SearchByTripletIPU.h"

SearchByTripletIPU::SearchByTripletIPU(
        poplar::Device device,
        std::vector<Event> events,
        bool is_model) :
        device(std::move(device)),
        events(events),
        graph(this->device),
        is_model(is_model)
{
    // Store the number of tiles.
    this->num_tiles = this->device.getTarget().getNumIPUs()
                    * this->device.getTarget().getTilesPerIPU();

    // Fix the number of threads to 3. We currently can't use all 6 due to
    // memory restrictions.
    this->num_threads = 3 * this->num_tiles;

    // Add codelets.
    this->graph.addCodelets({"src/SearchByTripletCodelet.cpp"}, "-O3 -I src");

    // Setup the graph program.
    this->setupGraphProgram();
}

std::tuple<std::vector<unsigned>, std::vector<Track>,
           std::vector<unsigned>, std::vector<Tracklet>>
    SearchByTripletIPU::execute(bool warmup, bool profile)
{
    // Warn user if they are using an IPUModel device.
    if (this->is_model)
    {
        std::cout << "Running on an IPUModel device. Benchmarking is meaningless!" << std::endl;
    }

    // Initialise the Poplar engine and load the IPU device.

    auto optionFlags = poplar::OptionFlags{};
    if (profile)
    {
        // Taken from UoB-HPC IPU cookbook:
        // https://github.com/UoB-HPC/ipu-hpc-cookbook
        optionFlags = poplar::OptionFlags
        {
            {"target.saveArchive",                "archive.a"},
            {"debug.instrument",                  "true"},
            {"debug.instrumentCompute",           "true"},
            {"debug.loweredVarDumpFile",          "vars.capnp"},
            {"debug.instrumentControlFlow",       "true"},
            {"debug.computeInstrumentationLevel", "tile"},
            {"debug.outputAllSymbols",            "true"},
            {"autoReport.all",                    "true"},
            {"autoReport.outputSerializedGraph",  "true"},
            {"debug.retainDebugInformation",      "true"}
        };
    }

    // Create buffers for IPU-to-host streams.
    std::vector<unsigned> num_tracks(this->num_threads);
    std::vector<unsigned> num_three_hit_tracks(this->num_threads);

    std::vector<Track> tracks(this->num_threads*constants::max_tracks);
    std::vector<Tracklet> three_hit_tracks(this->num_threads*constants::max_tracks);

    const unsigned track_size = this->num_threads* constants::max_tracks * (constants::max_track_size + 1);
    const unsigned three_hit_track_size = this->num_threads* constants::max_tracks * 3;

    // Create the engine and load the IPU device.
    poplar::Engine engine(this->graph, this->programs, optionFlags);
    engine.load(this->device);

    // Connect the data streams.
    engine.connectStream("module_pairs_write", this->module_pair_buffer.data());
    engine.connectStream("hits_write", this->hits_buffer.data());
    engine.connectStream("phi_write", this->phi_buffer.data());
    engine.connectStream("tracks_read", tracks.data());
    engine.connectStream("three_hit_tracks_read", three_hit_tracks.data());
    engine.connectStream("num_tracks_read", num_tracks.data());
    engine.connectStream("num_three_hit_tracks_read", num_tracks.data());

    // Perform a warmup run.
    if (warmup)
    {
        engine.run(Program::COPY_TO_IPU);
        engine.run(Program::RUN_ALGORITHM);
        engine.run(Program::COPY_FROM_IPU);
    }

    // Initialise the total run time in seconds.
    double total_secs = 0;

    // Record start time.
    auto start = std::chrono::high_resolution_clock::now();

    // Run the host-to-IPU data stream sequence.
    engine.run(Program::COPY_TO_IPU);

    // Record end time.
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

    // Calculate run time per event in seconds.
    auto secs = std::chrono::duration<double>(elapsed).count();

    // Update the total run time.
    total_secs += secs;

    // Report timing.
    std::cout << "Results...\n  Copy to IPU took: " << std::fixed << secs/1000 << " ms\n";

    // Reset start time.
    start = std::chrono::high_resolution_clock::now();

    // Run the main algorithm.
    engine.run(Program::RUN_ALGORITHM);

    // Record end time and work out run time.
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;

    // Work out run time in seconds.
    secs = std::chrono::duration<double>(elapsed).count();

    // Report timing.
    std::cout << "  Algorithm execution took: " << std::fixed << secs/1000 << " ms\n";
    auto events_per_sec = this->num_threads / secs;
    std::cout << "  Raw events per second: "<< std::fixed << events_per_sec << "\n";

    // Update the total run time.
    total_secs += secs;

    // Reset start time.
    start = std::chrono::high_resolution_clock::now();

    // Run the IPU-to-host data stream sequence.
    engine.run(Program::COPY_FROM_IPU);

    // Record end time and work out run time.
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;

    // Work out run time in seconds.
    secs = std::chrono::duration<double>(elapsed).count();

    // Update the total run time.
    total_secs += secs;

    // Report timing.
    std::cout << "  Copy from IPU took: " << std::fixed << secs/1000 << " ms\n";

    // Overall throughput.
    events_per_sec = this->num_threads / total_secs;
    std::cout << "  Events per second: " << std::fixed << events_per_sec << "\n";

    // Write profiling information to file.
    if (profile)
    {
        std::ofstream profile;
        profile.open("profile.txt");
        engine.printProfileSummary(profile, {{"showExecutionSteps", "true"}});
        profile.close();
    }

    return std::make_tuple(num_tracks, tracks,
                           num_three_hit_tracks, three_hit_tracks);
}

void SearchByTripletIPU::setupGraphProgram()
{
    // Store the number of events.
    auto num_events = this->events.size();

    // Store the number of module pairs. (Assume this is the same for all events.)
    auto num_module_pairs = this->events[0].getModulePairs().size();

    // Process one event per thread.
    if (num_events > this->num_threads)
    {
        throw std::runtime_error("Number of events must not exceed number of threads: "
            + std::to_string(num_threads));
    }

    // Size of buffer entries. (per tile)
    auto module_pair_size = 4 * constants::num_module_pairs;
    auto hits_size = 4 * constants::max_hits;
    auto phi_size =  constants::max_hits;
    auto candidates_size = constants::max_seeding_candidates;
    auto tracklets_size = constants::max_tracks_to_follow * 3;
    auto tracks_size = constants::max_tracks * (constants::max_track_size + 1);
    auto three_hit_tracks_size = constants::max_tracks * 3;
    auto track_mask_size = constants::max_tracks_to_follow;

    // Resize buffers.
    this->module_pair_buffer.resize(this->num_threads*module_pair_size);
    this->hits_buffer.resize(this->num_threads*hits_size);
    this->phi_buffer.resize(this->num_threads*phi_size);

    // Create the graph variables.

    poplar::Tensor module_pairs
        = this->graph.addVariable(poplar::FLOAT,
                                  {this->num_threads * module_pair_size},
                                  "module_pairs");
    poplar::Tensor hits
        = this->graph.addVariable(poplar::FLOAT,
                                  {this->num_threads * hits_size},
                                  "hits");
    poplar::Tensor phi
        = this->graph.addVariable(poplar::SHORT,
                                  {this->num_threads * phi_size},
                                  "phi");
    poplar::Tensor candidates
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_threads * candidates_size},
                                  "candidates");
    poplar::Tensor used_hits
        = this->graph.addVariable(poplar::BOOL,
                                  {this->num_threads * phi_size},
                                  "used_hits");
    poplar::Tensor tracklets
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_threads * tracklets_size},
                                  "tracklets");
    poplar::Tensor tracks
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_threads * tracks_size},
                                  "tracks");
    poplar::Tensor three_hit_tracks
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_threads * three_hit_tracks_size},
                                  "three_hit_tracks");
    poplar::Tensor num_tracks
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_threads},
                                  "num_tracks");
    poplar::Tensor num_three_hit_tracks
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_threads},
                                  "num_three_hit_tracks");
    poplar::Tensor track_mask
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_threads * track_mask_size},
                                  "track_mask");

    // Create a compute set.
    poplar::ComputeSet computeSet = this->graph.addComputeSet("computeSet");

    // Store the number of context workers (hardware threads) per tile.
    // We currently only use 3 of the 6 available threads due to memory
    // restrictions.
    const auto num_workers = 3;

    // Loop over the tiles, processing one event per thread.
    for (int i=0; i<this->num_threads; ++i)
    {
        // Work out the event index, modulo the number of threads.
        const auto event_idx = i%num_events;

        // Work out the tile index.
        const auto tile = int(i/num_workers);

        // Populate buffers and slice across threads.
        auto offset = i*module_pair_size;
        int idx = 0;
        for (const auto &module_pair : this->events[event_idx].getModulePairs())
        {
            module_pair_buffer[offset + 4*idx    ] = module_pair.hit_start;
            module_pair_buffer[offset + 4*idx + 1] = module_pair.num_hits;
            module_pair_buffer[offset + 4*idx + 2] = module_pair.z0;
            module_pair_buffer[offset + 4*idx + 3] = module_pair.z1;
            idx++;
        }
        this->graph.setTileMapping(module_pairs.slice(i*module_pair_size, (i+1)*module_pair_size), tile);

        offset = i*hits_size;
        idx = 0;
        for (const auto &hit : this->events[event_idx].getHits())
        {
            hits_buffer[offset + 4*idx    ] = hit.x;
            hits_buffer[offset + 4*idx + 1] = hit.y;
            hits_buffer[offset + 4*idx + 2] = hit.z;
            hits_buffer[offset + 4*idx + 3] = hit.module_idx;
            idx++;
        }
        this->graph.setTileMapping(hits.slice(i*hits_size, (i+1)*hits_size), tile);

        offset = i*phi_size;
        idx = 0;
        for (const auto &hit : this->events[event_idx].getHits())
        {
            phi_buffer[offset + idx] = hit.phi;
            idx++;
        }
        this->graph.setTileMapping(phi.slice(i*phi_size, (i+1)*phi_size), tile);
        this->graph.setTileMapping(candidates.slice(i*candidates_size, (i+1)*candidates_size), tile);
        this->graph.setTileMapping(used_hits.slice(i*phi_size, (i+1)*phi_size), tile);
        this->graph.setTileMapping(tracklets.slice(i*tracklets_size, (i+1)*tracklets_size), tile);
        this->graph.setTileMapping(tracks.slice(i*tracks_size, (i+1)*tracks_size), tile);
        this->graph.setTileMapping(three_hit_tracks.slice(i*three_hit_tracks_size, (i+1)*three_hit_tracks_size), tile);
        this->graph.setTileMapping(num_tracks.slice(i, i+1), tile);
        this->graph.setTileMapping(num_three_hit_tracks.slice(i, i+1), tile);
        this->graph.setTileMapping(track_mask.slice(i*track_mask_size, (i+1)*track_mask_size), tile);

        // Add a vertex to the compute set.
        poplar::VertexRef vtx = this->graph.addVertex(computeSet, "SearchByTriplet");

        // Connect variables to vertex inputs and outputs, sub-slicing tensors
        // over workers.
        this->graph.connect(
            vtx["module_pairs"],
            module_pairs.slice(i*module_pair_size, (i+1)*module_pair_size));
        this->graph.connect(
            vtx["hits"],
            hits.slice(i*hits_size, (i+1)*hits_size));
        this->graph.connect(
            vtx["used_hits"],
            used_hits.slice(i*phi_size, (i+1)*phi_size));
        this->graph.connect(
            vtx["phi"],
            phi.slice(i*phi_size, (i+1)*phi_size));
        this->graph.connect(
            vtx["candidate_hits"],
            candidates.slice(i*candidates_size, (i+1)*candidates_size));
        this->graph.connect(
            vtx["tracklets"],
            tracklets.slice(i*tracklets_size, (i+1)*tracklets_size));
        this->graph.connect(
            vtx["tracks"],
            tracks.slice(i*tracks_size, (i+1)*tracks_size));
        this->graph.connect(
            vtx["three_hit_tracks"],
            three_hit_tracks.slice(i*three_hit_tracks_size, (i+1)*three_hit_tracks_size));
        this->graph.connect(
            vtx["num_tracks"],
            num_tracks[i]);
        this->graph.connect(
            vtx["num_three_hit_tracks"],
            num_three_hit_tracks[i]);
        this->graph.connect(
            vtx["track_mask"],
            track_mask.slice(i*track_mask_size, (i+1)*track_mask_size));

        // Map the vertex to a tile.
        this->graph.setTileMapping(vtx, tile);

        // Add a performance estimate if this is an IPUModel device.
        if (this->is_model)
        {
            this->graph.setPerfEstimate(vtx, 100);
        }
    }

    // FIFO options. Allowing splitting of streams.
    auto fifo_opts = poplar::OptionFlags({{"splitLimit", "52428800"}});

    // Create host-to-IPU data streams and associated copy programs.

    auto module_pairs_write = this->graph.addHostToDeviceFIFO(
            "module_pairs_write",
            poplar::FLOAT,
            this->module_pair_buffer.size(),
            poplar::ReplicatedStreamMode::REPLICATE,
            fifo_opts);
    auto copy_module_pairs = poplar::program::Copy(module_pairs_write, module_pairs);

    auto hits_write = this->graph.addHostToDeviceFIFO(
            "hits_write",
            poplar::FLOAT,
            this->hits_buffer.size(),
            poplar::ReplicatedStreamMode::REPLICATE,
            fifo_opts);
    auto copy_hits = poplar::program::Copy(hits_write, hits);

    auto phi_write = this->graph.addHostToDeviceFIFO(
            "phi_write",
            poplar::SHORT,
            this->phi_buffer.size(),
            poplar::ReplicatedStreamMode::REPLICATE,
            fifo_opts);
    auto copy_phi = poplar::program::Copy(phi_write, phi);

    // Create a program sequence to copy data to the IPU.
    auto copy_to_ipu = poplar::program::Sequence
    {
        copy_module_pairs,
        copy_hits,
        copy_phi
    };

    // Add the host-to-IPU data transfer program.
    this->programs.push_back(copy_to_ipu);

    // Create IPU-to-host data streams and associated copy programs.

    const unsigned track_size = this->num_threads* constants::max_tracks * (constants::max_track_size + 1);
    const unsigned three_hit_track_size = this->num_threads* constants::max_tracks * 3;

    auto tracks_read = this->graph.addDeviceToHostFIFO(
            "tracks_read",
            poplar::UNSIGNED_INT,
            track_size,
            fifo_opts);
    auto copy_tracks = poplar::program::Copy(tracks, tracks_read);

    auto three_hit_tracks_read = this->graph.addDeviceToHostFIFO(
            "three_hit_tracks_read",
            poplar::UNSIGNED_INT,
            three_hit_track_size,
            fifo_opts);
    auto copy_three_hit_tracks = poplar::program::Copy(three_hit_tracks, three_hit_tracks_read);

    auto num_tracks_read = this->graph.addDeviceToHostFIFO(
            "num_tracks_read",
            poplar::UNSIGNED_INT,
            this->num_threads,
            fifo_opts);
    auto copy_num_tracks =
        poplar::program::Copy(num_tracks, num_tracks_read);

    auto num_three_hit_tracks_read = this->graph.addDeviceToHostFIFO(
            "num_three_hit_tracks_read",
            poplar::UNSIGNED_INT,
            this->num_threads,
            fifo_opts);
    auto copy_num_three_hit_tracks =
        poplar::program::Copy(num_three_hit_tracks, num_three_hit_tracks_read);

    // Add the compute set to the programs.
    this->programs.push_back(poplar::program::Execute(computeSet));

    // Create a program sequence to copy data from the IPU.
    auto copy_from_ipu = poplar::program::Sequence
    {
        copy_tracks,
        copy_three_hit_tracks,
        copy_num_tracks,
        copy_num_three_hit_tracks
    };

    // Add the IPU-to-host data transfer program.
    this->programs.push_back(copy_from_ipu);
}
