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

#include <popops/ElementWise.hpp>

#include "EventReader.h"
#include "SearchByTripletIPU.h"

// Set the number of workers per tile.
constexpr unsigned num_workers = 1;

SearchByTripletIPU::SearchByTripletIPU(
        poplar::Device device,
        std::vector<Event> events,
        const unsigned num_batches,
        bool ping_pong) :
        device(std::move(device)),
        events(events),
        graph(this->device),
        num_batches(num_batches),
        ping_pong(ping_pong)
{
    // Store the number of tiles (per IPU).
    this->num_tiles = this->device.getTarget().getTilesPerIPU();

    // Fix the number of threads to num_workers per tile. We currently can't
    // use all 6 due to memory restrictions.
    this->num_threads = num_workers*this->num_tiles;

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

    // Loop over threads to fill buffers.
    for (int i=0; i<this->num_threads; ++i)
    {
        // Work out the event index, modulo the number of threads.
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

    // Create buffers for IPU-to-host streams.
    std::vector<unsigned> num_tracks0(this->num_threads);
    std::vector<unsigned> num_tracks1(this->num_threads);
    std::vector<unsigned> num_three_hit_tracks0(this->num_threads);
    std::vector<unsigned> num_three_hit_tracks1(this->num_threads);

    std::vector<Track> tracks0(this->num_threads*constants::max_tracks);
    std::vector<Track> tracks1(this->num_threads*constants::max_tracks);
    std::vector<Tracklet> three_hit_tracks0(this->num_threads*constants::max_tracks);
    std::vector<Tracklet> three_hit_tracks1(this->num_threads*constants::max_tracks);

    // Create the engine and load the IPU device.
    poplar::Engine engine(this->graph, this->program, optionFlags);
    engine.load(this->device);

    // Copy data to the remote buffers.

    for (unsigned i=0; i<this->num_batches; ++i)
    {
        // Module pair records.
        engine.copyToRemoteBuffer(module_pair_buffer.data(), "module_pairs_rb0", i);
        engine.copyToRemoteBuffer(module_pair_buffer.data(), "module_pairs_rb1", i);
        engine.copyToRemoteBuffer(module_pair_buffer.data(), "module_pairs_rb2", i);
        engine.copyToRemoteBuffer(module_pair_buffer.data(), "module_pairs_rb3", i);

        // Raw hit data.
        engine.copyToRemoteBuffer(hits_buffer.data(), "hits_rb0", i);
        engine.copyToRemoteBuffer(hits_buffer.data(), "hits_rb1", i);
        engine.copyToRemoteBuffer(hits_buffer.data(), "hits_rb2", i);
        engine.copyToRemoteBuffer(hits_buffer.data(), "hits_rb3", i);

        // Raw phi data.
        engine.copyToRemoteBuffer(phi_buffer.data(), "phi_rb0", i);
        engine.copyToRemoteBuffer(phi_buffer.data(), "phi_rb1", i);
        engine.copyToRemoteBuffer(phi_buffer.data(), "phi_rb2", i);
        engine.copyToRemoteBuffer(phi_buffer.data(), "phi_rb3", i);
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

    auto events_per_sec = (4*this->num_threads*this->num_batches) / secs;
    std::cout << "Events per second: "<< std::fixed << events_per_sec << "\n";

    // Write profiling information to file.
    if (profile)
    {
        std::ofstream profile;
        profile.open("profile.txt");
        engine.printProfileSummary(profile, {{"showExecutionSteps", "true"}});
        profile.close();
    }

    std::cout << "Validating output...\n";

    // Copy data back from the remote buffers on IPUs 0 and 1 to validate that
    // the data is the same.

    // Number of tracks.
    engine.copyFromRemoteBuffer("num_tracks_rb0", num_tracks0.data(), 0);
    engine.copyFromRemoteBuffer("num_tracks_rb3", num_tracks1.data(), this->num_batches-1);
    for (unsigned i=0; i<num_tracks0.size(); ++i)
    {
        assert(num_tracks0[i] == num_tracks1[i]);
    }

    // Number of three-hit tracks.
    engine.copyFromRemoteBuffer("num_three_hit_tracks_rb0", num_three_hit_tracks0.data(), 0);
    engine.copyFromRemoteBuffer("num_three_hit_tracks_rb3", num_three_hit_tracks1.data(), this->num_batches-1);
    for (unsigned i=0; i<num_three_hit_tracks0.size(); ++i)
    {
        assert(num_three_hit_tracks0[i] == num_three_hit_tracks1[i]);
    }

    // Tracks.
    engine.copyFromRemoteBuffer("tracks_rb0", tracks0.data(), 0);
    engine.copyFromRemoteBuffer("tracks_rb3", tracks1.data(), this->num_batches-1);
    for (unsigned i=0; i<tracks0.size(); ++i)
    {
        assert(tracks0[i].num_hits == tracks1[i].num_hits);
        for (unsigned j=0; j<tracks0[i].num_hits; ++j)
        {
            assert(tracks0[i].hits[j] == tracks1[i].hits[j]);
        }
    }

    // Three-hit tracks.
    engine.copyFromRemoteBuffer("three_hit_tracks_rb0", three_hit_tracks0.data(), 0);
    engine.copyFromRemoteBuffer("three_hit_tracks_rb3", three_hit_tracks1.data(), this->num_batches-1);
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
    this->module_pair_buffer.resize(this->num_threads*module_pair_size);
    this->hits_buffer.resize(this->num_threads*hits_size);
    this->phi_buffer.resize(this->num_threads*phi_size);

    auto One = graph.addVariable(poplar::INT, {}, "1");
    this->graph.setTileMapping(One, 0);
    this->graph.setInitialValue(One, 1);

    // Set up indices for the repeat index of the remote buffers.
    auto rb_index0 = this->graph.addVariable(poplar::INT, {}, "rb_index0");
    auto rb_index1 = this->graph.addVariable(poplar::INT, {}, "rb_index1");
    auto rb_index2 = this->graph.addVariable(poplar::INT, {}, "rb_index2");
    auto rb_index3 = this->graph.addVariable(poplar::INT, {}, "rb_index3");
    this->graph.setTileMapping(rb_index0,                 0);
    this->graph.setTileMapping(rb_index1,   this->num_tiles);
    this->graph.setTileMapping(rb_index2, 2*this->num_tiles);
    this->graph.setTileMapping(rb_index3, 3*this->num_tiles);
    this->graph.setInitialValue(rb_index0, 0);
    this->graph.setInitialValue(rb_index1, 0);
    this->graph.setInitialValue(rb_index2, 0);
    this->graph.setInitialValue(rb_index3, 0);

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

    // Create the graph variables. For simplicity, we'll create duplicates
    // for each IPU and map them using the appropriate tile offset.

    // Module pair records.
    auto module_pairs0 = this->graph.addVariable(
        poplar::FLOAT, {this->num_threads*module_pair_size}, "module_pairs0");
    mapLinearlyOnOneIpu(module_pairs0, 0, this->device, this->graph);
    auto module_pairs_rb0 = this->graph.addRemoteBuffer(
        "module_pairs_rb0", poplar::FLOAT, this->num_threads*module_pair_size, this->num_batches);
    auto module_pairs1 = this->graph.addVariable(
        poplar::FLOAT, {this->num_threads*module_pair_size}, "module_pairs1");
    mapLinearlyOnOneIpu(module_pairs1, 1, this->device, this->graph);
    auto module_pairs_rb1 = this->graph.addRemoteBuffer(
        "module_pairs_rb1", poplar::FLOAT, this->num_threads*module_pair_size, this->num_batches);
    auto module_pairs2 = this->graph.addVariable(
        poplar::FLOAT, {this->num_threads*module_pair_size}, "module_pairs2");
    mapLinearlyOnOneIpu(module_pairs2, 2, this->device, this->graph);
    auto module_pairs_rb2 = this->graph.addRemoteBuffer(
        "module_pairs_rb2", poplar::FLOAT, this->num_threads*module_pair_size, this->num_batches);
    auto module_pairs3 = this->graph.addVariable(
        poplar::FLOAT, {this->num_threads*module_pair_size}, "module_pairs3");
    mapLinearlyOnOneIpu(module_pairs3, 3, this->device, this->graph);
    auto module_pairs_rb3 = this->graph.addRemoteBuffer(
        "module_pairs_rb3", poplar::FLOAT, this->num_threads*module_pair_size, this->num_batches);

    // Raw hits.
    auto hits0 = this->graph.addVariable(
        poplar::FLOAT, {this->num_threads*hits_size}, "hits0");
    mapLinearlyOnOneIpu(hits0, 0, this->device, this->graph);
    auto hits_rb0 = this->graph.addRemoteBuffer(
        "hits_rb0", poplar::FLOAT, this->num_threads*hits_size, this->num_batches);
    auto hits1 = this->graph.addVariable(
        poplar::FLOAT, {this->num_threads*hits_size}, "hits1");
    mapLinearlyOnOneIpu(hits1, 1, this->device, this->graph);
    auto hits_rb1 = this->graph.addRemoteBuffer(
        "hits_rb1", poplar::FLOAT, this->num_threads*hits_size, this->num_batches);
    auto hits2 = this->graph.addVariable(
        poplar::FLOAT, {this->num_threads*hits_size}, "hits2");
    mapLinearlyOnOneIpu(hits2, 2, this->device, this->graph);
    auto hits_rb2 = this->graph.addRemoteBuffer(
        "hits_rb2", poplar::FLOAT, this->num_threads*hits_size, this->num_batches);
    auto hits3 = this->graph.addVariable(
        poplar::FLOAT, {this->num_threads*hits_size}, "hits3");
    mapLinearlyOnOneIpu(hits3, 3, this->device, this->graph);
    auto hits_rb3 = this->graph.addRemoteBuffer(
        "hits_rb3", poplar::FLOAT, this->num_threads*hits_size, this->num_batches);

    // Phi values.
    auto phi0 = this->graph.addVariable(
        poplar::SHORT, {this->num_threads*phi_size}, "phi0");
    mapLinearlyOnOneIpu(phi0, 0, this->device, this->graph);
    auto phi_rb0 = this->graph.addRemoteBuffer(
        "phi_rb0", poplar::SHORT, this->num_threads*phi_size, this->num_batches);
    auto phi1 = this->graph.addVariable(
        poplar::SHORT, {this->num_threads*phi_size}, "phi1");
    mapLinearlyOnOneIpu(phi1, 1, this->device, this->graph);
    auto phi_rb1 = this->graph.addRemoteBuffer(
        "phi_rb1", poplar::SHORT, this->num_threads*phi_size, this->num_batches);
    auto phi2 = this->graph.addVariable(
        poplar::SHORT, {this->num_threads*phi_size}, "phi2");
    mapLinearlyOnOneIpu(phi2, 2, this->device, this->graph);
    auto phi_rb2 = this->graph.addRemoteBuffer(
        "phi_rb2", poplar::SHORT, this->num_threads*phi_size, this->num_batches);
    auto phi3 = this->graph.addVariable(
        poplar::SHORT, {this->num_threads*phi_size}, "phi3");
    mapLinearlyOnOneIpu(phi3, 3, this->device, this->graph);
    auto phi_rb3 = this->graph.addRemoteBuffer(
        "phi_rb3", poplar::SHORT, this->num_threads*phi_size, this->num_batches);

    // Candidate hits for a track.
    auto candidates0 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*candidates_size}, "candidates0");
    mapLinearlyOnOneIpu(candidates0, 0, this->device, this->graph);
    auto candidates1 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*candidates_size}, "candidates1");
    mapLinearlyOnOneIpu(candidates1, 1, this->device, this->graph);
    auto candidates2 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*candidates_size}, "candidates2");
    mapLinearlyOnOneIpu(candidates2, 2, this->device, this->graph);
    auto candidates3 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*candidates_size}, "candidates3");
    mapLinearlyOnOneIpu(candidates3, 3, this->device, this->graph);

    // Hits that have been assigned to a track.
    auto used_hits0 = this->graph.addVariable(
        poplar::BOOL, {this->num_threads*phi_size}, "used_hits0");
    mapLinearlyOnOneIpu(used_hits0, 0, this->device, this->graph);
    auto used_hits1 = this->graph.addVariable(
        poplar::BOOL, {this->num_threads*phi_size}, "used_hits1");
    mapLinearlyOnOneIpu(used_hits1, 1, this->device, this->graph);
    auto used_hits2 = this->graph.addVariable(
        poplar::BOOL, {this->num_threads*phi_size}, "used_hits2");
    mapLinearlyOnOneIpu(used_hits2, 2, this->device, this->graph);
    auto used_hits3 = this->graph.addVariable(
        poplar::BOOL, {this->num_threads*phi_size}, "used_hits3");
    mapLinearlyOnOneIpu(used_hits3, 3, this->device, this->graph);

    // Seed tracks.
    auto tracklets0 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*tracklets_size}, "tracklets0");
    mapLinearlyOnOneIpu(tracklets0, 0, this->device, this->graph);
    auto tracklets1 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*tracklets_size}, "tracklets1");
    mapLinearlyOnOneIpu(tracklets1, 1, this->device, this->graph);
    auto tracklets2 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*tracklets_size}, "tracklets2");
    mapLinearlyOnOneIpu(tracklets2, 2, this->device, this->graph);
    auto tracklets3 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*tracklets_size}, "tracklets3");
    mapLinearlyOnOneIpu(tracklets3, 3, this->device, this->graph);

    // Assigned tracks.
    auto tracks0 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*tracks_size}, "tracks0");
    mapLinearlyOnOneIpu(tracks0, 0, this->device, this->graph);
    auto tracks_rb0 = this->graph.addRemoteBuffer(
        "tracks_rb0", poplar::UNSIGNED_INT, this->num_threads*tracks_size, this->num_batches);
    auto tracks1 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*tracks_size}, "tracks1");
    mapLinearlyOnOneIpu(tracks1, 1, this->device, this->graph);
    auto tracks_rb1 = this->graph.addRemoteBuffer(
        "tracks_rb1", poplar::UNSIGNED_INT, this->num_threads*tracks_size, this->num_batches);
    auto tracks2 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*tracks_size}, "tracks2");
    mapLinearlyOnOneIpu(tracks2, 2, this->device, this->graph);
    auto tracks_rb2 = this->graph.addRemoteBuffer(
        "tracks_rb2", poplar::UNSIGNED_INT, this->num_threads*tracks_size, this->num_batches);
    auto tracks3 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*tracks_size}, "tracks3");
    mapLinearlyOnOneIpu(tracks3, 3, this->device, this->graph);
    auto tracks_rb3 = this->graph.addRemoteBuffer(
        "tracks_rb3", poplar::UNSIGNED_INT, this->num_threads*tracks_size, this->num_batches);

    // Assigned three-hit tracks.
    auto three_hit_tracks0 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*three_hit_tracks_size}, "three_hit_tracks0");
    mapLinearlyOnOneIpu(three_hit_tracks0, 0, this->device, this->graph);
    auto three_hit_tracks_rb0 = this->graph.addRemoteBuffer(
        "three_hit_tracks_rb0", poplar::UNSIGNED_INT, this->num_threads*three_hit_tracks_size, this->num_batches);
    auto three_hit_tracks1 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*three_hit_tracks_size}, "three_hit_tracks1");
    mapLinearlyOnOneIpu(three_hit_tracks1, 1, this->device, this->graph);
    auto three_hit_tracks_rb1 = this->graph.addRemoteBuffer(
        "three_hit_tracks_rb1", poplar::UNSIGNED_INT, this->num_threads*three_hit_tracks_size, this->num_batches);
    auto three_hit_tracks2 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*three_hit_tracks_size}, "three_hit_tracks2");
    mapLinearlyOnOneIpu(three_hit_tracks2, 2, this->device, this->graph);
    auto three_hit_tracks_rb2 = this->graph.addRemoteBuffer(
        "three_hit_tracks_rb2", poplar::UNSIGNED_INT, this->num_threads*three_hit_tracks_size, this->num_batches);
    auto three_hit_tracks3 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*three_hit_tracks_size}, "three_hit_tracks3");
    mapLinearlyOnOneIpu(three_hit_tracks3, 3, this->device, this->graph);
    auto three_hit_tracks_rb3 = this->graph.addRemoteBuffer(
        "three_hit_tracks_rb3", poplar::UNSIGNED_INT, this->num_threads*three_hit_tracks_size, this->num_batches);

    // Number of assigned tracks.
    auto num_tracks0 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads}, "num_tracks0");
    mapLinearlyOnOneIpu(num_tracks0, 0, this->device, this->graph);
    auto num_tracks_rb0 = this->graph.addRemoteBuffer(
        "num_tracks_rb0", poplar::UNSIGNED_INT, this->num_threads, this->num_batches);
    auto num_tracks1 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads}, "num_tracks1");
    mapLinearlyOnOneIpu(num_tracks1, 1, this->device, this->graph);
    auto num_tracks_rb1 = this->graph.addRemoteBuffer(
        "num_tracks_rb1", poplar::UNSIGNED_INT, this->num_threads, this->num_batches);
    auto num_tracks2 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads}, "num_tracks2");
    mapLinearlyOnOneIpu(num_tracks2, 2, this->device, this->graph);
    auto num_tracks_rb2 = this->graph.addRemoteBuffer(
        "num_tracks_rb2", poplar::UNSIGNED_INT, this->num_threads, this->num_batches);
    auto num_tracks3 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads}, "num_tracks3");
    mapLinearlyOnOneIpu(num_tracks3, 3, this->device, this->graph);
    auto num_tracks_rb3 = this->graph.addRemoteBuffer(
        "num_tracks_rb3", poplar::UNSIGNED_INT, this->num_threads, this->num_batches);

    // Number of assigned three-hit tracks.
    auto num_three_hit_tracks0 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads}, "num_three_hit_tracks0");
    mapLinearlyOnOneIpu(num_three_hit_tracks0, 0, this->device, this->graph);
    auto num_three_hit_tracks_rb0 = this->graph.addRemoteBuffer(
        "num_three_hit_tracks_rb0", poplar::UNSIGNED_INT, this->num_threads, this->num_batches);
    auto num_three_hit_tracks1 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads}, "num_three_hit_tracks1");
    mapLinearlyOnOneIpu(num_three_hit_tracks1, 1, this->device, this->graph);
    auto num_three_hit_tracks_rb1 = this->graph.addRemoteBuffer(
        "num_three_hit_tracks_rb1", poplar::UNSIGNED_INT, this->num_threads, this->num_batches);
    auto num_three_hit_tracks2 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads}, "num_three_hit_tracks2");
    mapLinearlyOnOneIpu(num_three_hit_tracks2, 2, this->device, this->graph);
    auto num_three_hit_tracks_rb2 = this->graph.addRemoteBuffer(
        "num_three_hit_tracks_rb2", poplar::UNSIGNED_INT, this->num_threads, this->num_batches);
    auto num_three_hit_tracks3 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads}, "num_three_hit_tracks3");
    mapLinearlyOnOneIpu(num_three_hit_tracks3, 3, this->device, this->graph);
    auto num_three_hit_tracks_rb3 = this->graph.addRemoteBuffer(
        "num_three_hit_tracks_rb3", poplar::UNSIGNED_INT, this->num_threads, this->num_batches);

    // Track identity mask.
    auto track_mask0 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*track_mask_size}, "track_mask0");
    mapLinearlyOnOneIpu(track_mask0, 0, this->device, this->graph);
    auto track_mask1 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*track_mask_size}, "track_mask1");
    mapLinearlyOnOneIpu(track_mask1, 1, this->device, this->graph);
    auto track_mask2 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*track_mask_size}, "track_mask2");
    mapLinearlyOnOneIpu(track_mask2, 2, this->device, this->graph);
    auto track_mask3 = this->graph.addVariable(
        poplar::UNSIGNED_INT, {this->num_threads*track_mask_size}, "track_mask3");
    mapLinearlyOnOneIpu(track_mask3, 3, this->device, this->graph);

    // Lambda to create the algorithm program for each IPU.
    auto createAlgorithmProgram = [&](
            poplar::Tensor &module_pairs,
            poplar::Tensor &hits,
            poplar::Tensor &phi,
            poplar::Tensor &candidates,
            poplar::Tensor &used_hits,
            poplar::Tensor &tracklets,
            poplar::Tensor &tracks,
            poplar::Tensor &three_hit_tracks,
            poplar::Tensor &num_tracks,
            poplar::Tensor &num_three_hit_tracks,
            poplar::Tensor &track_mask,
            const int ipu_num)
    {
        // Create a compute set.
        auto computeSet = this->graph.addComputeSet("computeSet");

        // Add a vertex for each thread that will be used on the tile
        // and connect tensor slices to vertex inputs and outputs.

        // Work out the starting tile index.
        const auto start_tile = ipu_num*this->num_tiles;

        // Loop over all threads on the tile.
        for (unsigned i=0; i<this->num_threads; ++i)
        {
            // Work out the tile index.
            const auto tile = start_tile + int(i/num_workers);

            // Add a vertex to the compute set.
            auto vtx = this->graph.addVertex(computeSet, "SearchByTriplet");

            // Connect variables to vertex inputs and outputs, sub-slicing
            // tensors over workers.
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

        return poplar::program::Execute(computeSet);
    };

    // Create the algorithm execution program for each IPU.

    // Run on IPU 0.
    const auto execute_algorithm_ipu0 = poplar::program::Sequence
    {
        createAlgorithmProgram(
            module_pairs0,
            hits0,
            phi0,
            candidates0,
            used_hits0,
            tracklets0,
            tracks0,
            three_hit_tracks0,
            num_tracks0,
            num_three_hit_tracks0,
            track_mask0,
            0
        )
    };

    // Run on IPU 1.
    const auto execute_algorithm_ipu1 = poplar::program::Sequence
    {
        createAlgorithmProgram(
            module_pairs1,
            hits1,
            phi1,
            candidates1,
            used_hits1,
            tracklets1,
            tracks1,
            three_hit_tracks1,
            num_tracks1,
            num_three_hit_tracks1,
            track_mask1,
            1
        )
    };

    // Run on IPU 2.
    const auto execute_algorithm_ipu2 = poplar::program::Sequence
    {
        createAlgorithmProgram(
            module_pairs2,
            hits2,
            phi2,
            candidates2,
            used_hits2,
            tracklets2,
            tracks2,
            three_hit_tracks2,
            num_tracks2,
            num_three_hit_tracks2,
            track_mask2,
            2
        )
    };

    // Run on IPU 3.
    const auto execute_algorithm_ipu3 = poplar::program::Sequence
    {
        createAlgorithmProgram(
            module_pairs3,
            hits3,
            phi3,
            candidates3,
            used_hits3,
            tracklets3,
            tracks3,
            three_hit_tracks3,
            num_tracks3,
            num_three_hit_tracks3,
            track_mask3,
            3
        )
    };

    // Create programs to copy data to/from the IPUs using the remote buffers.

    // Module pair records.
    const auto copy_from_module_pairs_rb0 =
        poplar::program::Copy(module_pairs_rb0, module_pairs0, rb_index0);
    const auto copy_from_module_pairs_rb1 =
        poplar::program::Copy(module_pairs_rb1, module_pairs1, rb_index1);
    const auto copy_from_module_pairs_rb2 =
        poplar::program::Copy(module_pairs_rb2, module_pairs2, rb_index2);
    const auto copy_from_module_pairs_rb3 =
        poplar::program::Copy(module_pairs_rb3, module_pairs3, rb_index3);

    // Raw hits.
    const auto copy_from_hits_rb0 =
        poplar::program::Copy(hits_rb0, hits0, rb_index0);
    const auto copy_from_hits_rb1 =
        poplar::program::Copy(hits_rb1, hits1, rb_index1);
    const auto copy_from_hits_rb2 =
        poplar::program::Copy(hits_rb2, hits2, rb_index2);
    const auto copy_from_hits_rb3 =
        poplar::program::Copy(hits_rb3, hits3, rb_index3);

    // Phi values.
    const auto copy_from_phi_rb0 =
        poplar::program::Copy(phi_rb0, phi0, rb_index0);
    const auto copy_from_phi_rb1 =
        poplar::program::Copy(phi_rb1, phi1, rb_index1);
    const auto copy_from_phi_rb2 =
        poplar::program::Copy(phi_rb2, phi2, rb_index2);
    const auto copy_from_phi_rb3 =
        poplar::program::Copy(phi_rb3, phi3, rb_index3);

    // Assigned tracks.
    const auto copy_to_tracks_rb0 =
        poplar::program::Copy(tracks0, tracks_rb0, rb_index0);
    const auto copy_to_tracks_rb1 =
        poplar::program::Copy(tracks1, tracks_rb1, rb_index1);
    const auto copy_to_tracks_rb2 =
        poplar::program::Copy(tracks2, tracks_rb2, rb_index2);
    const auto copy_to_tracks_rb3 =
        poplar::program::Copy(tracks3, tracks_rb3, rb_index3);

    // Assigned three-hit tracks.
    const auto copy_to_three_hit_tracks_rb0 =
        poplar::program::Copy(three_hit_tracks0, three_hit_tracks_rb0, rb_index0);
    const auto copy_to_three_hit_tracks_rb1 =
        poplar::program::Copy(three_hit_tracks1, three_hit_tracks_rb1, rb_index1);
    const auto copy_to_three_hit_tracks_rb2 =
        poplar::program::Copy(three_hit_tracks2, three_hit_tracks_rb2, rb_index2);
    const auto copy_to_three_hit_tracks_rb3 =
        poplar::program::Copy(three_hit_tracks3, three_hit_tracks_rb3, rb_index3);

    // Number of assigned tracks.
    const auto copy_to_num_tracks_rb0 =
        poplar::program::Copy(num_tracks0, num_tracks_rb0, rb_index0);
    const auto copy_to_num_tracks_rb1 =
        poplar::program::Copy(num_tracks1, num_tracks_rb1, rb_index1);
    const auto copy_to_num_tracks_rb2 =
        poplar::program::Copy(num_tracks2, num_tracks_rb2, rb_index2);
    const auto copy_to_num_tracks_rb3 =
        poplar::program::Copy(num_tracks3, num_tracks_rb3, rb_index3);

    // Number of assigned three-hit tracks.
    const auto copy_to_num_three_hit_tracks_rb0 =
        poplar::program::Copy(num_three_hit_tracks0, num_three_hit_tracks_rb0, rb_index0);
    const auto copy_to_num_three_hit_tracks_rb1 =
        poplar::program::Copy(num_three_hit_tracks1, num_three_hit_tracks_rb1, rb_index1);
    const auto copy_to_num_three_hit_tracks_rb2 =
        poplar::program::Copy(num_three_hit_tracks2, num_three_hit_tracks_rb2, rb_index2);
    const auto copy_to_num_three_hit_tracks_rb3 =
        poplar::program::Copy(num_three_hit_tracks3, num_three_hit_tracks_rb3, rb_index3);

    // Consolidate copy programs for each IPU.

    // To IPU 0.
    const auto copy_to_ipu0 = poplar::program::Sequence
    {
        copy_from_module_pairs_rb0,
        copy_from_hits_rb0,
        copy_from_phi_rb0
    };

    // To IPU 1.
    const auto copy_to_ipu1 = poplar::program::Sequence
    {
        copy_from_module_pairs_rb1,
        copy_from_hits_rb1,
        copy_from_phi_rb1
    };

    // To IPU 3.
    const auto copy_to_ipu2 = poplar::program::Sequence
    {
        copy_from_module_pairs_rb2,
        copy_from_hits_rb2,
        copy_from_phi_rb2
    };

    // To IPU 3.
    const auto copy_to_ipu3 = poplar::program::Sequence
    {
        copy_from_module_pairs_rb3,
        copy_from_hits_rb3,
        copy_from_phi_rb3
    };

    // From IPU 0.
    const auto copy_from_ipu0 = poplar::program::Sequence
    {
        copy_to_tracks_rb0,
        copy_to_three_hit_tracks_rb0,
        copy_to_num_tracks_rb0,
        copy_to_num_three_hit_tracks_rb0
    };

    // From IPU 1.
    const auto copy_from_ipu1 = poplar::program::Sequence
    {
        copy_to_tracks_rb1,
        copy_to_three_hit_tracks_rb1,
        copy_to_num_tracks_rb1,
        copy_to_num_three_hit_tracks_rb1
    };

    // From IPU 2.
    const auto copy_from_ipu2 = poplar::program::Sequence
    {
        copy_to_tracks_rb2,
        copy_to_three_hit_tracks_rb2,
        copy_to_num_tracks_rb2,
        copy_to_num_three_hit_tracks_rb2
    };

    // From IPU 3.
    const auto copy_from_ipu3 = poplar::program::Sequence
    {
        copy_to_tracks_rb3,
        copy_to_three_hit_tracks_rb3,
        copy_to_num_tracks_rb3,
        copy_to_num_three_hit_tracks_rb3
    };

    // Lambda function to increment the remote buffer indices.
    const auto increment = [&](poplar::Tensor &t)
    {
        poplar::program::Sequence s;
        popops::addInPlace(graph, t, One, s, "t++");
        return s;
    };

    if (this->ping_pong)
    {
        // Run in ping-pong fashion, i.e. alternating compute and data transfer
        // between pairs of IPU devices. This means that IPUs 0 and 1 compute
        // while IPUs 2 and 3 load data from the exchange. Next, IPUs 2 and 3
        // compute while IPUs 0 and 1 copy results back to the exchange.
        this->program = poplar::program::Sequence
        {
            poplar::program::Sequence
            {
                copy_to_ipu0,
                copy_to_ipu1
            },
            poplar::program::Sequence
            {
                execute_algorithm_ipu0,
                execute_algorithm_ipu1,
                copy_to_ipu2,
                copy_to_ipu3
            },
            poplar::program::Sequence
            {
                copy_from_ipu0,
                copy_from_ipu1,
                execute_algorithm_ipu2,
                execute_algorithm_ipu3
            },
            poplar::program::Repeat(
                this->num_batches-1,
                poplar::program::Sequence
                {
                    poplar::program::Sequence
                    {
                        increment(rb_index0),
                        increment(rb_index1),
                    },
                    poplar::program::Sequence
                    {
                        copy_to_ipu0,
                        copy_to_ipu1,
                        copy_from_ipu2,
                        copy_from_ipu3
                    },
                    poplar::program::Sequence
                    {
                        increment(rb_index2),
                        increment(rb_index3),
                    },
                    poplar::program::Sequence
                    {
                        execute_algorithm_ipu0,
                        execute_algorithm_ipu1,
                        copy_to_ipu2,
                        copy_to_ipu3
                    },
                    poplar::program::Sequence
                    {
                        copy_from_ipu0,
                        copy_from_ipu1,
                        execute_algorithm_ipu2,
                        execute_algorithm_ipu3,
                    },
                }
            ),
            poplar::program::Sequence
            {
                copy_from_ipu2,
                copy_from_ipu3
            }
        };
    }
    else
    {
        // Run in a regular fashion, alternativing compute and data transfer
        // for all of the IPUs, i.e. first copy data from the exchange to each
        // IPU, then compute on the IPUs, and finally copy data from the IPUs
        // back to the exchange.
        this->program = poplar::program::Sequence
        {
            poplar::program::Repeat(
                this->num_batches,
                poplar::program::Sequence
                {
                    poplar::program::Sequence
                    {
                        copy_to_ipu0,
                        copy_to_ipu1,
                        copy_to_ipu2,
                        copy_to_ipu3
                    },
                    poplar::program::Sequence
                    {
                        execute_algorithm_ipu0,
                        execute_algorithm_ipu1,
                        execute_algorithm_ipu2,
                        execute_algorithm_ipu3
                    },
                    poplar::program::Sequence
                    {
                        copy_from_ipu0,
                        copy_from_ipu1,
                        copy_from_ipu2,
                        copy_from_ipu3
                    },
                    poplar::program::Sequence
                    {
                        increment(rb_index0),
                        increment(rb_index1),
                        increment(rb_index2),
                        increment(rb_index3),
                    },
                }
            ),
        };
    }
}
