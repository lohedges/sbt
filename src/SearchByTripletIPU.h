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

#ifndef _SEARCHBYTRIPLETIPU_H
#define _SEARCHBYTRIPLETIPU_H

#include<tuple>
#include<vector>

#include <poplar/DeviceManager.hpp>
#include <poplar/Engine.hpp>
#include <poplar/Graph.hpp>

#include "Definitions.h"

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

        \parm is_model
            Whether the IPU device is an IPUModel.
     */
     SearchByTripletIPU(poplar::Device device,
                        std::vector<Event> events,
                        bool is_model=false);

    //! Execute the Kalman filter.
    /*! \param events_per_sec
            The number of events processed per second.

        \param warmup
            Whether to perform a "warmup" run. This is useful when benchmarking.

        \param profile
            Whether to write profiling information to file.

        \return num_tracks, tracks, num_three_hit_tracks, tracklets
            The reconstructed tracks and three-hit tracks.
     */
    std::tuple<std::vector<unsigned>, std::vector<Track>,
               std::vector<unsigned>, std::vector<Tracklet>>
        execute(
            double&,
            bool warmup=false,
            bool profile=false);

private:
    /// The vector of events.
    std::vector<Event> events;

    /// The Poplar device.
    poplar::Device device;

    /// The Poplar graph.
    poplar::Graph graph;

    /// The Poplar program sequence.
    poplar::program::Sequence prog;

    /// The number of IPU tiles.
    unsigned num_tiles;

    /// Whether the device is an IPUModel.
    bool is_model;

    /// Set up the graph program.
    void setupGraphProgram();

    /// Buffer for storing flattened module pair data.
    std::vector<float> module_pair_buffer;

    /// Buffer for storing flattened hit data.
    std::vector<float> hits_buffer;

    /// Buffer for storing flattened phi data.
    std::vector<int16_t> phi_buffer;
};

#endif	/* _SEARCHBYTRIPLETIPU_H */
