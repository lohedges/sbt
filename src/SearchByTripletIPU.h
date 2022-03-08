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

#ifndef _SEARCHBYTRIPLETIPU_H
#define _SEARCHBYTRIPLETIPU_H

#include <vector>

#include <poplar/DeviceManager.hpp>
#include <poplar/Engine.hpp>
#include <poplar/Graph.hpp>

#include "Definitions.h"

enum Program
{
    COPY_TO_IPU,
    RUN_ALGORITHM,
    COPY_FROM_IPU
};

struct Tracklet
{
    unsigned hits[3];
};

struct Track
{
    unsigned num_hits;
    unsigned hits[constants::max_track_size];
};

class SearchByTripletIPU
{
public:
    //! Constructor.
    /*! \param device
            The Poplar device.

        \param events
            A vector of events.

        \param num_batches
            The number of batches to process.
     */
     SearchByTripletIPU(poplar::Device device,
                        std::vector<Event> events,
                        const unsigned num_batches=50);

    //! Execute the Search by Triplet algorithm and report timing statistics.
    /*! \param warmup
            Whether to perform a "warmup" run. This can be useful when
            benchmarking.

        \param profile
            Whether to write profiling information to file.
     */
    void execute(bool warmup=false, bool profile=false);

private:
    /// The vector of events.
    std::vector<Event> events;

    /// The Poplar device.
    poplar::Device device;

    /// The Poplar graph.
    poplar::Graph graph;

    /// The Poplar program.
    poplar::program::Program program;

	/// The Poplar program sequence.
    std::vector<poplar::program::Program> programs;

    /// The number of batches to process.
    unsigned num_batches;

    /// The number of IPUs.
    unsigned num_ipus;

    /// The number of tiles per IPU.
    unsigned num_tiles;

    /// The number of threads per tile.
    unsigned num_threads;

    /// Create the graph program using data streams.
    void createGraphProgramDataStreams();

    /// Create the graph program using remote buffers streams.
    void createGraphProgramRemoteBuffers();

    //! Execute the Search by Triplet algorithm using data streams and report
    //  timing statistics.
    /*! \param warmup
            Whether to perform a "warmup" run. This can be useful when
            benchmarking.

        \param profile
            Whether to write profiling information to file.
     */
    void executeDataStreams(bool warmup=false, bool profile=false);

    //! Execute the Search by Triplet algorithm using remote buffers and report
    //  timing statistics.
    /*! \param warmup
            Whether to perform a "warmup" run. This can be useful when
            benchmarking.

        \param profile
            Whether to write profiling information to file.
     */
    void executeRemoteBuffers(bool warmup=false, bool profile=false);

    /// Buffer for storing flattened module pair data.
    std::vector<float> module_pair_buffer;

    /// Buffer for storing flattened hit data.
    std::vector<float> hits_buffer;

    /// Buffer for storing flattened phi data.
    std::vector<int16_t> phi_buffer;
};

#endif	/* _SEARCHBYTRIPLETIPU_H */
