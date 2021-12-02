/*
  Copyright (c) 2021 Lester Hedges <lester.hedges@gmail.com>

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
    this->num_tiles = this->device.getTarget().getTilesPerIPU();

    // Add codelets.
    this->graph.addCodelets({"src/SearchByTripletCodelet.cpp"}, "-O3 -I src");

    // Setup the graph program.
    this->setupGraphProgram();
}

std::tuple<unsigned*, Track*, unsigned*, Tracklet*> SearchByTripletIPU::execute(
    double &events_per_sec,
    bool warmup,
    bool profile)
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

    poplar::Engine engine(this->graph, this->prog, optionFlags);
    engine.load(this->device);

    engine.writeTensor("module_pairs_write",
      this->module_pair_buffer.data(),
      this->module_pair_buffer.data() + this->module_pair_buffer.size());

    engine.writeTensor("hits_write",
      this->hits_buffer.data(),
      this->hits_buffer.data() + this->hits_buffer.size());

    engine.writeTensor("phi_write",
      this->phi_buffer.data(),
      this->phi_buffer.data() + this->phi_buffer.size());

    // Perform a warmup run.
    if (warmup)
        engine.run(0);

    // Record start time.
    auto start = std::chrono::high_resolution_clock::now();

    // Run the program.
    engine.run(0);

    // Record end time.
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;

    // Calculate run time per event in seconds.
    auto secs = std::chrono::duration<double>(elapsed).count();
    events_per_sec = this->num_tiles / secs;

    // Write profiling information to file.
    if (profile)
    {
        std::ofstream profile;
        profile.open("profile.txt");
        engine.printProfileSummary(profile, {{"showExecutionSteps", "true"}});
        profile.close();
    }

    unsigned *num_tracks = new unsigned[this->num_tiles];
    unsigned *num_three_hit_tracks = new unsigned[this->num_tiles];

    Track *tracks = new Track[this->num_tiles*constants::max_tracks];
    Tracklet *three_hit_tracks = new Tracklet[this->num_tiles*constants::max_tracks];

    unsigned track_size = this->num_tiles* constants::max_tracks * (constants::max_track_size + 1);
    unsigned three_hit_track_size = this->num_tiles* constants::max_tracks * 3;

    engine.readTensor("num_tracks_read", num_tracks, num_tracks + this->num_tiles);
    engine.readTensor("num_three_hit_tracks_read", num_three_hit_tracks, num_three_hit_tracks + this->num_tiles);
    engine.readTensor("tracks_read", tracks, tracks + track_size);
    engine.readTensor("three_hit_tracks_read", three_hit_tracks, three_hit_tracks + three_hit_track_size);

    return std::make_tuple(num_tracks, tracks,
                           num_three_hit_tracks, three_hit_tracks);
}

void SearchByTripletIPU::setupGraphProgram()
{
    // Store the number of events.
    auto num_events = this->events.size();

    // Store the number of module pairs. (Assume this is the same for all events.)
    auto num_module_pairs = this->events[0].getModulePairs().size();

    // Process one event per tile if there are less events than IPUs.
    if (num_events > this->num_tiles)
    {
        throw std::runtime_error("Number of events must not exceed number of tiles: "
            + std::to_string(num_tiles));
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
    this->module_pair_buffer.resize(this->num_tiles*module_pair_size);
    this->hits_buffer.resize(this->num_tiles*hits_size);
    this->phi_buffer.resize(this->num_tiles*phi_size);

    // Create the graph variables.

    poplar::Tensor module_pairs
        = this->graph.addVariable(poplar::FLOAT,
                                  {this->num_tiles * module_pair_size},
                                  "module_pairs");
    poplar::Tensor hits
        = this->graph.addVariable(poplar::FLOAT,
                                  {this->num_tiles * hits_size},
                                  "hits");
    poplar::Tensor phi
        = this->graph.addVariable(poplar::SHORT,
                                  {this->num_tiles * phi_size},
                                  "phi");
    poplar::Tensor candidates
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_tiles * candidates_size},
                                  "candidates");
    poplar::Tensor used_hits
        = this->graph.addVariable(poplar::BOOL,
                                  {this->num_tiles * phi_size},
                                  "used_hits");
    poplar::Tensor tracklets
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_tiles * tracklets_size},
                                  "tracklets");
    poplar::Tensor tracks
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_tiles * tracks_size},
                                  "tracks");
    poplar::Tensor three_hit_tracks
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_tiles * three_hit_tracks_size},
                                  "three_hit_tracks");
    poplar::Tensor num_tracks
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_tiles},
                                  "num_tracks");
    poplar::Tensor num_three_hit_tracks
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_tiles},
                                  "num_three_hit_tracks");
    poplar::Tensor track_mask
        = this->graph.addVariable(poplar::UNSIGNED_INT,
                                  {this->num_tiles * track_mask_size},
                                  "track_mask");

    // Create a compute set.
    poplar::ComputeSet computeSet = this->graph.addComputeSet("computeSet");

    // Loop over the tiles, processing one event per tile.
    for (int i=0; i<this->num_tiles; ++i)
    {
        // Work out the event index, modulo the number of tiles.
        const auto event_idx = i%num_events;

        // Populate buffers and slice across tiles.
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
        this->graph.setTileMapping(module_pairs.slice(i*module_pair_size, (i+1)*module_pair_size), i);

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
        this->graph.setTileMapping(hits.slice(i*hits_size, (i+1)*hits_size), i);

        offset = i*phi_size;
        idx = 0;
        for (const auto &hit : this->events[event_idx].getHits())
        {
            phi_buffer[offset + idx] = hit.phi;
            idx++;
        }
        this->graph.setTileMapping(phi.slice(i*phi_size, (i+1)*phi_size), i);
        this->graph.setTileMapping(candidates.slice(i*candidates_size, (i+1)*candidates_size), i);
        this->graph.setTileMapping(used_hits.slice(i*phi_size, (i+1)*phi_size), i);
        this->graph.setTileMapping(tracklets.slice(i*tracklets_size, (i+1)*tracklets_size), i);
        this->graph.setTileMapping(tracks.slice(i*tracks_size, (i+1)*tracks_size), i);
        this->graph.setTileMapping(three_hit_tracks.slice(i*three_hit_tracks_size, (i+1)*three_hit_tracks_size), i);
        this->graph.setTileMapping(num_tracks.slice(i, i+1), i);
        this->graph.setTileMapping(num_three_hit_tracks.slice(i, i+1), i);
        this->graph.setTileMapping(track_mask.slice(i*track_mask_size, (i+1)*track_mask_size), i);

        // Create the vertex identification string.
        std::string vertex_id = std::to_string(i);

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
        this->graph.setTileMapping(vtx, i);

        // Add a performance estimate if this is an IPUModel device.
        if (this->is_model)
        {
            this->graph.setPerfEstimate(vtx, 100);
        }
    }

    // Create host read/write handles for buffers.
    this->graph.createHostWrite("module_pairs_write", module_pairs);
    this->graph.createHostWrite("hits_write", hits);
    this->graph.createHostWrite("phi_write", phi);
    this->graph.createHostWrite("used_hits_write", used_hits);
    this->graph.createHostRead("tracks_read", tracks);
    this->graph.createHostRead("three_hit_tracks_read", three_hit_tracks);
    this->graph.createHostRead("num_tracks_read", num_tracks);
    this->graph.createHostRead("num_three_hit_tracks_read", num_three_hit_tracks);

    // Add the compute set to the program.
    this->prog.add(poplar::program::Execute(computeSet));
}
