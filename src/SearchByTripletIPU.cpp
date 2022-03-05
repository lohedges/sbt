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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <popops/ElementWise.hpp>

#include "EventReader.h"
#include "SearchByTripletIPU.h"

constexpr unsigned num_threads = 1;

SearchByTripletIPU::SearchByTripletIPU(
        poplar::Device device,
        std::vector<Event> events,
        const unsigned num_batches) :
        device(std::move(device)),
        events(events),
        graph(this->device),
        num_batches(num_batches)
{
    // Store the number of IPUs on the target device.
    this->num_ipus = this->device.getTarget().getNumIPUs();

    // Store the number of tiles (per IPU).
    this->num_tiles = this->device.getTarget().getTilesPerIPU();

    // Add codelets.
    this->graph.addCodelets({"src/SearchByTripletCodelet.cpp"}, "-O3 -I src");

    // Setup the graph program.
    this->setupGraphProgram();
}

void  SearchByTripletIPU::execute(bool warmup, bool profile)
{
    // Store the number of events.
    const auto num_events = this->events.size();

    // Size of buffer entries. (per tile)
    const auto module_pair_size = 4*constants::num_module_pairs;
    const auto hits_size = 4*constants::max_hits;
    const auto phi_size = constants::max_hits;

    // Loop over tiles to fill buffers.
    for (int i=0; i<this->num_tiles; ++i)
    {
        // Work out the event index.
        const auto event_idx = i%num_events;

        // Populate buffers.
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

        offset = i*phi_size;
        idx = 0;
        for (const auto &hit : this->events[event_idx].getHits())
        {
            phi_buffer[offset + idx] = hit.phi;
            idx++;
        }
    }

    // Initialise the Poplar engine and load the IPU device.

    auto optionFlags = poplar::OptionFlags{};
    if (profile)
    {
        // Taken from UoB-HPC IPU cookbook:
        // https://github.com/UoB-HPC/ipu-hpc-cookbook
        optionFlags = poplar::OptionFlags
        {
            {"target.extendedMemory",             "true"},
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
    else
    {
        optionFlags = poplar::OptionFlags
        {
            {"target.extendedMemory", "true"},
        };
    }

    // Create the engine and load the IPU device.
    poplar::Engine engine(this->graph, this->program, optionFlags);
    engine.load(this->device);

    // Copy data to the remote buffers.
    std::cout << "Copying data to remote buffers...\n";

    // Loop over batches.
    for (unsigned i=0; i<this->num_batches; ++i)
    {
        // Loop over tiles.
        for (unsigned j=0; j<this->num_ipus; ++j)
        {
            // Loop over threads.
            for (unsigned k=0; k<num_threads; ++k)
            {
                // Module pair records.
                std::string name = "module_pairs_rb_" + std::to_string(j)
                                 + std::to_string(k);
                engine.copyToRemoteBuffer(module_pair_buffer.data(), name, i);

                // Raw hit data.
                name = "hits_rb_" + std::to_string(j) + std::to_string(k);
                engine.copyToRemoteBuffer(hits_buffer.data(), name, i);

                // Raw phi data.
                name = "phi_rb_" + std::to_string(j) + std::to_string(k);
                engine.copyToRemoteBuffer(phi_buffer.data(), name, i);
            }
        }
    }

    std::cout << "Running benchmarks...\n";

    // Record start time.
    auto start = std::chrono::high_resolution_clock::now();

    // Run program.
    engine.run(0);

    // Record end time.
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

    // Calculate run time per event in seconds.
    auto secs = std::chrono::duration<double>(elapsed).count();

    auto events_per_sec = (this->num_ipus*num_threads*this->num_tiles*this->num_batches) / secs;
    std::cout << "Events per second: "<< std::fixed << events_per_sec << "\n";

    // Write profiling information to file.
    if (profile)
    {
        std::ofstream profile;
        profile.open("profile.txt");
        engine.printProfileSummary(profile, {{"showExecutionSteps", "true"}});
        profile.close();
    }

    // Copy data back from the remote buffers to validate the the data is
    // the same for different IPUs, threads, and batches. We validate data
    // on IPU 0, thread 0, batch 0 against that of num_ipus-1, num_thread-1,
    // num_batches-1.
    std::cout << "Validating output...\n";

    // Create buffers for IPU-to-host streams.
    std::vector<unsigned> num_tracks0(this->num_tiles);
    std::vector<unsigned> num_tracks1(this->num_tiles);
    std::vector<unsigned> num_three_hit_tracks0(this->num_tiles);
    std::vector<unsigned> num_three_hit_tracks1(this->num_tiles);

    std::vector<Track> tracks0(this->num_tiles*constants::max_tracks);
    std::vector<Track> tracks1(this->num_tiles*constants::max_tracks);
    std::vector<Tracklet> three_hit_tracks0(this->num_tiles*constants::max_tracks);
    std::vector<Tracklet> three_hit_tracks1(this->num_tiles*constants::max_tracks);

    // Number of tracks.
    std::string name = "num_tracks_rb_" + std::to_string(this->num_ipus-1)
                                        + std::to_string(num_threads-1);
    engine.copyFromRemoteBuffer("num_tracks_rb_00", num_tracks0.data(), 0);
    engine.copyFromRemoteBuffer(name, num_tracks1.data(), this->num_batches-1);
    for (unsigned i=0; i<num_tracks0.size(); ++i)
    {
        assert(num_tracks0[i] == num_tracks1[i]);
    }

    // Number of three-hit tracks.
    name = "num_three_hit_tracks_rb_" + std::to_string(this->num_ipus-1)
                                      + std::to_string(num_threads-1);
    engine.copyFromRemoteBuffer("num_three_hit_tracks_rb_00", num_three_hit_tracks0.data(), 0);
    engine.copyFromRemoteBuffer(name, num_three_hit_tracks1.data(), this->num_batches-1);
    for (unsigned i=0; i<num_three_hit_tracks0.size(); ++i)
    {
        assert(num_three_hit_tracks0[i] == num_three_hit_tracks1[i]);
    }

    // Tracks.
    name = "tracks_rb_" + std::to_string(this->num_ipus-1)
                        + std::to_string(num_threads-1);
    engine.copyFromRemoteBuffer("tracks_rb_00", tracks0.data(), 0);
    engine.copyFromRemoteBuffer(name, tracks1.data(), this->num_batches-1);
    for (unsigned i=0; i<tracks0.size(); ++i)
    {
        assert(tracks0[i].num_hits == tracks1[i].num_hits);
        for (unsigned j=0; j<tracks0[i].num_hits; ++j)
        {
            assert(tracks0[i].hits[j] == tracks1[i].hits[j]);
        }
    }

    // Three-hit tracks.
    name = "three_hit_tracks_rb_" + std::to_string(this->num_ipus-1)
                                  + std::to_string(num_threads-1);
    engine.copyFromRemoteBuffer("three_hit_tracks_rb_00", three_hit_tracks0.data(), 0);
    engine.copyFromRemoteBuffer(name, three_hit_tracks1.data(), this->num_batches-1);
    for (unsigned i=0; i<three_hit_tracks0.size(); ++i)
    {
        assert(three_hit_tracks0[i].hits[0] == three_hit_tracks1[i].hits[0]);
        assert(three_hit_tracks0[i].hits[1] == three_hit_tracks1[i].hits[1]);
        assert(three_hit_tracks0[i].hits[2] == three_hit_tracks1[i].hits[2]);
    }
}

void SearchByTripletIPU::setupGraphProgram()
{
    // We will use remote buffers to hold num_batches copies of the event data.
    // This will be processed by 2 IPUs at a time in a ping-pong fashion, i.e.
    // IPUs 0 and 1 will compute while IPUs 2 and 3 transfer data, and vice-versa.

    // Store the number of module pairs. (Assume this is the same for all events.)
    const auto num_module_pairs = this->events[0].getModulePairs().size();

    // Size of buffer entries. (per tile)
    const auto module_pair_size = 4*constants::num_module_pairs;
    const auto hits_size = 4*constants::max_hits;
    const auto phi_size = constants::max_hits;
    const auto candidates_size = constants::max_seeding_candidates;
    const auto tracklets_size = constants::max_tracks_to_follow * 3;
    const auto tracks_size = constants::max_tracks * (constants::max_track_size + 1);
    const auto three_hit_tracks_size = constants::max_tracks * 3;
    const auto track_mask_size = constants::max_tracks_to_follow;

    // Resize buffers.
    this->module_pair_buffer.resize(this->num_tiles*module_pair_size);
    this->hits_buffer.resize(this->num_tiles*hits_size);
    this->phi_buffer.resize(this->num_tiles*phi_size);

    auto One = graph.addVariable(poplar::INT, {}, "1");
    this->graph.setTileMapping(One, 0);
    this->graph.setInitialValue(One, 1);

    // Set up indices for the repeat index of the remote buffers.
    auto rb_index = this->graph.addVariable(poplar::INT, {}, "rb_index");
    this->graph.setTileMapping(rb_index, 0);

    // Taken from UoB-HPC IPU cookbook:
    // https://github.com/UoB-HPC/ipu-hpc-cookbook
    // This is required to ensure that tensors only exist on a single IPU,
    // which is not the case when using Poplar's mapLinearlyWithOffset.
    const auto mapLinearlyOnOneIpu = [&](
        poplar::Tensor &tensor,
        const int ipuNum,
        poplar::Device &device,
        poplar::Graph &graph)
    {
        auto totalElements = 1;
        for (auto dim: tensor.shape())
        {
            totalElements *= dim;
        }
        int numTilesPerIpu = device.getTarget().getNumTiles() / device.getTarget().getNumIPUs();
        auto itemsPerTile = totalElements / numTilesPerIpu;
        auto numTilesWithExtraItem = totalElements % numTilesPerIpu;
        auto from = 0;
        for (auto tileNum = 0; tileNum < numTilesPerIpu; tileNum++)
        {
            const auto itemsForThisTile = (tileNum < numTilesWithExtraItem) ? itemsPerTile + 1 : itemsPerTile;
            const auto to = from + itemsForThisTile;
            graph.setTileMapping(tensor.slice(from, to), tileNum + (ipuNum * numTilesPerIpu));
            from = to;
        }
    };

    // Lambda to create the programs for each IPU and thread.
    auto createPrograms = [&](
        const int ipu_num,
        const int thread_num,
        poplar::program::Sequence &compute,
        poplar::program::Sequence &copy_to_ipu,
        poplar::program::Sequence &copy_from_ipu)
    {
        // Module pair records.
        std::string name = "module_pairs_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto module_pairs = this->graph.addVariable(
            poplar::FLOAT, {this->num_tiles*module_pair_size}, name);
        mapLinearlyOnOneIpu(module_pairs, ipu_num, this->device, this->graph);
        name = "module_pairs_rb_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto module_pairs_rb = this->graph.addRemoteBuffer(
            name, poplar::FLOAT, this->num_tiles*module_pair_size, this->num_batches);

        // Raw hits.
        name = "hits_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto hits = this->graph.addVariable(
            poplar::FLOAT, {this->num_tiles*hits_size}, name);
        mapLinearlyOnOneIpu(hits, ipu_num, this->device, this->graph);
        name = "hits_rb_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto hits_rb = this->graph.addRemoteBuffer(
            name, poplar::FLOAT, this->num_tiles*hits_size, this->num_batches);

        // Phi values.
        name = "phi_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto phi = this->graph.addVariable(
            poplar::SHORT, {this->num_tiles*phi_size}, name);
        mapLinearlyOnOneIpu(phi, ipu_num, this->device, this->graph);
        name = "phi_rb_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto phi_rb = this->graph.addRemoteBuffer(
            name, poplar::SHORT, this->num_tiles*phi_size, this->num_batches);

        // Caldidate hits.
        name = "candidates_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto candidates = this->graph.addVariable(
            poplar::UNSIGNED_INT, {this->num_tiles*candidates_size}, name);
        mapLinearlyOnOneIpu(candidates, ipu_num, this->device, this->graph);

        // Hits that have been assigned to tracks.
        name = "used_hits_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto used_hits = this->graph.addVariable(
            poplar::BOOL, {this->num_tiles*phi_size}, name);
        mapLinearlyOnOneIpu(used_hits, ipu_num, this->device, this->graph);

        // Seed tracks.
        name = "tracklets_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto tracklets = this->graph.addVariable(
            poplar::UNSIGNED_INT, {this->num_tiles*tracklets_size}, name);
        mapLinearlyOnOneIpu(tracklets, ipu_num, this->device, this->graph);

        // Assigned tracks.
        name = "tracks_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto tracks = this->graph.addVariable(
            poplar::UNSIGNED_INT, {this->num_tiles*tracks_size}, name);
        mapLinearlyOnOneIpu(tracks, ipu_num, this->device, this->graph);
        name = "tracks_rb_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto tracks_rb = this->graph.addRemoteBuffer(
            name, poplar::UNSIGNED_INT, this->num_tiles*tracks_size, this->num_batches);

        // Assigned three-hit tracks.
        name = "three_hit_tracks_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto three_hit_tracks = this->graph.addVariable(
            poplar::UNSIGNED_INT, {this->num_tiles*three_hit_tracks_size}, name);
        mapLinearlyOnOneIpu(three_hit_tracks, ipu_num, this->device, this->graph);
        name = "three_hit_tracks_rb_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto three_hit_tracks_rb = this->graph.addRemoteBuffer(
            name, poplar::UNSIGNED_INT, this->num_tiles*three_hit_tracks_size, this->num_batches);

        // Number of tracks.
        name = "num_tracks_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto num_tracks = this->graph.addVariable(
            poplar::UNSIGNED_INT, {this->num_tiles}, name);
        mapLinearlyOnOneIpu(num_tracks, ipu_num, this->device, this->graph);
        name = "num_tracks_rb_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto num_tracks_rb = this->graph.addRemoteBuffer(
            name, poplar::UNSIGNED_INT, this->num_tiles, this->num_batches);

        // Number of three-hit tracks.
        name = "num_three_hit_tracks_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto num_three_hit_tracks = this->graph.addVariable(
            poplar::UNSIGNED_INT, {this->num_tiles}, name);
        mapLinearlyOnOneIpu(num_three_hit_tracks, ipu_num, this->device, this->graph);
        name = "num_three_hit_tracks_rb_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto num_three_hit_tracks_rb = this->graph.addRemoteBuffer(
            name, poplar::UNSIGNED_INT, this->num_tiles, this->num_batches);

        // Track identity mask.
        name = "track_mask_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto track_mask = this->graph.addVariable(
            poplar::UNSIGNED_INT, {this->num_tiles*track_mask_size}, name);
        mapLinearlyOnOneIpu(track_mask, ipu_num, this->device, this->graph);

        // Create a compute set.
        name = "compute_set_" + std::to_string(ipu_num) + std::to_string(thread_num);
        auto compute_set = this->graph.addComputeSet(name);

        // Work out the starting tile index.
        const auto start_tile = ipu_num*this->num_tiles;

        // Loop over all tiles and add a vertex for each, connecting inputs
        // and outputs to the appropriate slices of the tensors above.
        for (unsigned i=0; i<this->num_tiles; ++i)
        {
            // Work out the tile index.
            const auto tile = start_tile + i;

            // Add a vertex to the compute set.
            auto vtx = this->graph.addVertex(compute_set, "SearchByTriplet");

            // Connect variables to vertex inputs and outputs.
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

            // Map the vertex to the tile.
            this->graph.setTileMapping(vtx, tile);
        }

        // Add to the compute programs for this tile.
        compute.add(poplar::program::Execute(compute_set));

        // Add to the copy to IPU programs for this tile.
        copy_to_ipu.add(poplar::program::Copy(module_pairs_rb, module_pairs, rb_index));
        copy_to_ipu.add(poplar::program::Copy(hits_rb, hits, rb_index));
        copy_to_ipu.add(poplar::program::Copy(phi_rb, phi, rb_index));

        // Add to the copy from IPU programs for this tile.
        copy_from_ipu.add(poplar::program::Copy(tracks, tracks_rb, rb_index));
        copy_from_ipu.add(poplar::program::Copy(three_hit_tracks, three_hit_tracks_rb, rb_index));
        copy_from_ipu.add(poplar::program::Copy(num_tracks, num_tracks_rb, rb_index));
        copy_from_ipu.add(poplar::program::Copy(num_three_hit_tracks, num_three_hit_tracks_rb, rb_index));
    };

    // Compute program sequences for each stage, i.e. copy to IPU, compute, and
    // copy from IPU.
    poplar::program::Sequence compute;
    poplar::program::Sequence copy_to_ipu;
    poplar::program::Sequence copy_from_ipu;

    // Loop over all IPUs.
    for (unsigned i=0; i<this->num_ipus; ++i)
    {
        // Loop over all threads.
        for (unsigned j=0; j<num_threads; ++j)
        {
            // Add programs for each IPU and thread.
            createPrograms(i, j, compute, copy_to_ipu, copy_from_ipu);
        }
    }

    // Lambda function to increment the remote buffer index.
    const auto increment = [&](poplar::Tensor &t)
    {
        poplar::program::Sequence s;
        popops::addInPlace(graph, t, One, s, "t++");
        return s;
    };

    // Run alternating compute and data transferfor all of the IPUs, i.e.
    // first copy data from the exchange to each IPU, then compute on the IPUs,
    // and finally copy data from the IPUs back to the exchange.
    this->program = poplar::program::Sequence
    {
        poplar::program::Repeat(
            this->num_batches,
            poplar::program::Sequence
            {
                copy_to_ipu,
                compute,
                copy_from_ipu,
                increment(rb_index)
            }
        ),
    };
}
