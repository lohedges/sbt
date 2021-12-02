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

// Adapted from the LHCb Allen project: https://gitlab.cern.ch/lhcb/Allen

#include <tuple>

#include <poplar/Vertex.hpp>

#include "Utils.h"

using ModulePair = struct
{
    float hit_start;
    float num_hits;
    float z0;
    float z1;
};

using Tracklet = struct
{
    unsigned hits[3];
};

using Track = struct
{
    unsigned num_hits;
    unsigned hits[constants::max_track_size];
};

using Hit = struct
{
    float x;
    float y;
    float z;
    float module_idx;
};

using namespace poplar;

// Handy aliases for Poplar Input and InOut types. (These are rank-1 tensors.)
// https://docs.graphcore.ai/projects/assembly-programming/en/latest/vertex_vectors.html
// Use default alignment for types.
using InputShortTensor    = Input<Vector<int16_t>>;
using InputFloatTensor    = Input<Vector<float>>;
using InOutUnsignedTensor = InOut<Vector<unsigned>>;
using InOutBoolTensor     = InOut<Vector<bool>>;

// Helper function prototypes.

// Find track seeding candidates for a given module pair triplet.
void track_seeding(
    ModulePair *module_pairs,
    const unsigned mp_idx,
    Hit *hits,
    Tracklet *tracklets,
    unsigned &tracks_to_follow,
    int16_t *phi,
    bool *is_used,
    unsigned *track_mask,
    unsigned *initial_seeding_candidates);

// Forward the candidates to the next module pair triplet.
void track_forwarding(
    ModulePair *module_pairs,
    const unsigned mp_idx,
    Hit *hits,
    Tracklet *tracklets,
    Track *tracks,
    Tracklet *three_hit_tracks,
    unsigned &tracks_to_follow,
    unsigned &num_tracks,
    unsigned &num_three_hit_tracks,
    int16_t *phi,
    bool *is_used,
    unsigned *track_mask,
    const unsigned diff_ttf,
    const unsigned prev_ttf);

// Find first candidate in the forwards direction.
std::tuple<int16_t, int16_t> find_forward_candidate(
    ModulePair &module_pair,
    int16_t *phi,
    const Hit &h0,
    const float tx,
    const float ty,
    const float dz);

// Search by triplet vertex.
class SearchByTriplet : public Vertex
{
public:
    // Input fields. (constant)
    InputFloatTensor module_pairs;              // Module pair records.
    InputFloatTensor hits;                      // Packed hits for all module pairs. {x, y, z, ...}
    InputShortTensor phi;                       // Packed phi values for all hits.

    // InOut fields. (read/writeable)
    InOutBoolTensor     used_hits;              // Whether a hit has been assigned to a track.
    InOutUnsignedTensor candidate_hits;         // Candidate hits for a track.
    InOutUnsignedTensor tracklets;              // Three-hit tracks used for seeding.
    InOutUnsignedTensor tracks;                 // The assigned tracks.
    InOutUnsignedTensor three_hit_tracks;       // The assigned three-hit tracks.
    InOut<unsigned>     num_tracks;             // The number of tracks.
    InOut<unsigned>     num_three_hit_tracks;   // The number of three-hit tracks.
    InOutUnsignedTensor track_mask;             // The masked track numbers.

    // Overloaded compute function.
    bool compute() {
        // Re-cast attributes for convenience.

        // Cast the module_pairs tensor as an array of ModulePair structs.
        auto module_pairs_array =
            *static_cast<::ModulePair(*)[constants::max_track_size]>(static_cast<void*>(&module_pairs[0]));

        // Cast the hits tensor as an array of Hit structs.
        auto hit_array =
            *static_cast<::Hit(*)[constants::max_hits]>(static_cast<void*>(&hits[0]));

        // Cast the tracklets tensor as an array of Track structs.
        auto tracklet_array =
            *static_cast<::Tracklet(*)[constants::max_tracks_to_follow]>(static_cast<void*>(&tracklets[0]));

        // Cast the tracks tensor as an array of Track structs.
        auto track_array =
            *static_cast<::Track(*)[constants::max_tracks]>(static_cast<void*>(&tracks[0]));

        // Cast the three-hit tracks tensor as an array of Tracklet structs.
        auto three_hit_track_array =
            *static_cast<::Tracklet(*)[constants::max_tracks]>(static_cast<void*>(&three_hit_tracks[0]));

        // Cast the phi tensor as an array of int16 types.
        auto phi_array =
            *static_cast<int16_t(*)[phi.size()]>(static_cast<void*>(&phi[0]));

        // Cast the used_hits tensor as an array of bool types.
        auto used_hit_array =
            *static_cast<bool(*)[used_hits.size()]>(static_cast<void*>(&used_hits[0]));

        // Cast the track mask tensor as an array of unsigned ints.
        auto track_mask_array =
            *static_cast<unsigned(*)[track_mask.size()]>(static_cast<void*>(&track_mask[0]));

        // Cast the candidates tensor as an array of unsigned ints.
        auto candidate_array =
            *static_cast<unsigned(*)[candidate_hits.size()]>(static_cast<void*>(&candidate_hits[0]));

        // Initialise variables.
        unsigned tracks_to_follow = 0;
        unsigned last_ttf = 0;
        unsigned first_module_pair = constants::max_track_size - 2;

        // Zero the track counters.
        *num_tracks = 0;
        *num_three_hit_tracks = 0;

        // Set all hits as unused.
        memset(used_hit_array, false, used_hits.size() * sizeof(bool));

        // Perform intitial track seeding of first module pair triplet.
        // Start at the last three module pairs.
        track_seeding(
            module_pairs_array,
            first_module_pair,
            hit_array,
            tracklet_array,
            tracks_to_follow,
            phi_array,
            used_hit_array,
            track_mask_array,
            candidate_array);

        --first_module_pair;

        // Forwarding / seeding loop over remaining module pair triplets,
        // working backwards through the module pairs.
        while (first_module_pair > 0)
        {
            const auto prev_ttf = last_ttf;
            last_ttf = tracks_to_follow;
            const auto diff_ttf = last_ttf - prev_ttf;

            // Track forwarding.
            track_forwarding(
                module_pairs_array,
                first_module_pair,
                hit_array,
                tracklet_array,
                track_array,
                three_hit_track_array,
                tracks_to_follow,
                *num_tracks,
                *num_three_hit_tracks,
                phi_array,
                used_hit_array,
                track_mask_array,
                diff_ttf,
                prev_ttf);

            // Track seeding.
            track_seeding(
                module_pairs_array,
                first_module_pair,
                hit_array,
                tracklet_array,
                tracks_to_follow,
                phi_array,
                used_hit_array,
                track_mask_array,
                candidate_array);

            --first_module_pair;
        }

        const auto prev_ttf = last_ttf;
        last_ttf = tracks_to_follow;
        const auto diff_ttf = last_ttf - prev_ttf;

        // Process the last bunch of tracks.
        for (unsigned i=0; i<diff_ttf; ++i)
        {
            const auto full_track_number = track_mask[(prev_ttf + i) % constants::max_tracks_to_follow];
            const bool track_flag = (full_track_number & bits::seed) == bits::seed;

            // Here we are only interested in three-hit tracks.
            if (track_flag)
            {
                const auto track_number = full_track_number & bits::track_number;
                const Tracklet* tracklet = tracklet_array + track_number;

                three_hit_track_array[num_three_hit_tracks].hits[0] = tracklet->hits[0];
                three_hit_track_array[num_three_hit_tracks].hits[1] = tracklet->hits[1];
                three_hit_track_array[num_three_hit_tracks].hits[2] = tracklet->hits[2];

                num_three_hit_tracks++;
            }
        }

        return true;
    }
};

// Helper function definitions.

void track_seeding(
    ModulePair *module_pairs,
    const unsigned mp_idx,
    Hit *hits,
    Tracklet *tracklets,
    unsigned &tracks_to_follow,
    int16_t *phi,
    bool *is_used,
    unsigned *track_mask,
    unsigned *initial_seeding_candidates)
{
    // Store the hit offset for this module pair.
    auto offset = int(module_pairs[mp_idx].hit_start);

    // Loop over all hits for this module pair.
    for (unsigned i=0; i<module_pairs[mp_idx].num_hits; ++i)
    {
        // Get the index of the hit in the phi array.
        auto h1_index = i + offset;

        // Make sure the hit hasn't already been used.
        if (not is_used[h1_index])
        {
            // Initialise the outputs.
            uint16_t best_h0 = 0;
            uint16_t best_h2 = 0;
            float best_fit = constants::max_scatter;

            // Extract the hit.
            const auto h1 = hits[h1_index];

            // Get the phi value of the hit.
            const auto h1_phi = phi[h1_index];

            // Find the first candidate in the previous module.
            auto phi_index = binary_search_leftmost(
                phi + int(module_pairs[mp_idx+1].hit_start),
                module_pairs[mp_idx+1].num_hits,
                h1_phi);

            // Initialise the number of candidates that have been found.
            unsigned found_h0_candidates = 0;

            // Do a "pendulum search" to find candidates, consisting of iterating in the
            // following manner:
            // phi_index, phi_index + 1, phi_index - 1, phi_index + 2, ...
            int num_hits = int(module_pairs[mp_idx+1].num_hits);
            for (unsigned j=0;
                 j < num_hits && found_h0_candidates < constants::max_seeding_candidates;
                 ++j)
            {
                // By setting the sign to the oddity of j, the search behaviour
                // is achieved.
                const auto sign = j & 0x01;
                const int index_diff = sign ? j : -j;
                phi_index += index_diff;

                const auto index_in_bounds =
                    (phi_index < 0 ? phi_index +  num_hits :
                                    (phi_index >= num_hits ?
                                     phi_index -  num_hits :
                                     phi_index));

                const auto h0_index = int(module_pairs[mp_idx+1].hit_start) + index_in_bounds;

                // Only store the candidate if it hasn't been used.
                if (not is_used[h0_index])
                {
                    initial_seeding_candidates[found_h0_candidates++] = h0_index;
                }
            }

            // Now use the candidates found previously (initial_seeding_candidates) to find the best
            // triplet. Since data is sorted, search using a binary search.
            for (unsigned j=0; j<found_h0_candidates; ++j)
            {
                // Get the hit index and extract the hit.
                const auto h0_index = initial_seeding_candidates[j];
                const auto h0 = hits[h0_index];

                // Initialise extrapolation parameters.
                const auto td = 1.0f / (h1.z - h0.z);
                const auto txn = (h1.x - h0.x);
                const auto tyn = (h1.y - h0.y);
                const auto tx = txn * td;
                const auto ty = tyn * td;

                // Get the candidates by performing a binary search in expected phi.
                const auto candidate_h2 = find_forward_candidate(
                    module_pairs[mp_idx-1],
                    phi,
                    h0,
                    tx,
                    ty,
                    module_pairs[mp_idx-1].z0 - module_pairs[mp_idx+1].z0);

                // Extract the first candidate in the next module pair.
                // Since the buffer is circular, finding the container size means finding
                // the first element.
                const auto candidate_h2_index = std::get<0>(candidate_h2);
                const auto extrapolated_phi = std::get<1>(candidate_h2);

                const int num_hits = int(module_pairs[mp_idx-1].num_hits);
                for (unsigned k=0; k<num_hits; ++k)
                {
                    const auto index_in_bounds = (candidate_h2_index + k) % num_hits;
                    const auto h2_index = int(module_pairs[mp_idx-1].hit_start) + index_in_bounds;

                    // Check that the phi difference is within the tolerance using
                    // modulo arithmentic.
                    const int16_t phi_diff = phi[h2_index] - extrapolated_phi;
                    const int16_t abs_phi_diff = phi_diff < 0 ? -phi_diff : phi_diff;
                    if (abs_phi_diff > constants::phi_tolerance_i16)
                    {
                        break;
                    }

                    if (not is_used[h2_index])
                    {
                        const auto h2 = hits[h2_index];
                        const auto dz = h2.z - h0.z;
                        const auto predx = h0.x + tx * dz;
                        const auto predy = h0.y + ty * dz;
                        const auto dx = predx - h2.x;
                        const auto dy = predy - h2.y;

                        const auto scatter = dx*dx + dy*dy;

                        // Is this the best fit to date?
                        if (scatter < best_fit)
                        {
                            best_fit = scatter;
                            best_h0 = h0_index;
                            best_h2 = h2_index;
                        }
                    }
                }
            }

            if (best_fit < constants::max_scatter)
            {
                // Work out the track number, modulo the maximum number of tracks.
                const auto track_number = tracks_to_follow % constants::max_tracks_to_follow;

                // Store the tracklet.
                tracklets[track_number].hits[0] = best_h0;
                tracklets[track_number].hits[1] = h1_index;
                tracklets[track_number].hits[2] = best_h2;

                // Store the track mask.
                track_mask[track_number] = bits::seed | track_number;

                // Update the number of tracks to follow.
                tracks_to_follow++;
            }
        }
    }
}

void track_forwarding(
    ModulePair *module_pairs,
    const unsigned mp_idx,
    Hit *hits,
    Tracklet *tracklets,
    Track *tracks,
    Tracklet *three_hit_tracks,
    unsigned &tracks_to_follow,
    unsigned &num_tracks,
    unsigned &num_three_hit_tracks,
    int16_t *phi,
    bool *is_used,
    unsigned *track_mask,
    const unsigned diff_ttf,
    const unsigned prev_ttf)
{
    for (unsigned i=0; i<diff_ttf; ++i)
    {
        // Get the full masked track number.
        const auto full_track_number = track_mask[(prev_ttf + i) % constants::max_tracks_to_follow];

        // Extract information about the track.
        const bool track_flag = (full_track_number & bits::seed) == bits::seed;
        const auto skipped_modules = (full_track_number & bits::skipped_modules) >> bits::skipped_module_position;
        auto track_number = full_track_number & bits::track_number;

        Hit h0, h1;
        Track *track;
        Tracklet *tracklet;
        unsigned number_of_hits;

        if (track_flag)
        {
            tracklet = tracklets + track_number;
            number_of_hits = 3;

            h0 = hits[tracklet->hits[number_of_hits - 2]];
            h1 = hits[tracklet->hits[number_of_hits - 1]];
        }
        else
        {
            track = tracks + track_number;
            number_of_hits = track->num_hits;

            h0 = hits[track->hits[number_of_hits - 2]];
            h1 = hits[track->hits[number_of_hits - 1]];
        }

        // Get the z coordinate for h0 within the module pair.
        float z;
        if (int(h0.module_idx) % 2 == 0)
        {
            z = module_pairs[mp_idx-1].z0;
        }
        else
        {
            z = module_pairs[mp_idx-1].z1;
        }

        // Track forwarding over t, for all hits in the next module.
        // Line calculations.
        const auto td  = 1.0f / (h1.z - h0.z);
        const auto txn = (h1.x - h0.x);
        const auto tyn = (h1.y - h0.y);
        const auto tx  = txn * td;
        const auto ty  = tyn * td;

        // Find the best candidate.
        float best_fit = constants::max_scatter;
        int best_h2 = -1;

        // Get candidates by performing a binary search in expected phi.
        const auto candidate_h2 = find_forward_candidate(
                    module_pairs[mp_idx-1],
                    phi,
                    h0,
                    tx,
                    ty,
                    z - h0.z);

        // First candidate in the next module pair.
        const auto candidate_h2_index = std::get<0>(candidate_h2);
        const auto extrapolated_phi = std::get<1>(candidate_h2);

        int num_hits = int(module_pairs[mp_idx-1].num_hits);
        for (unsigned int j=0; j<num_hits; ++j)
        {
            const auto index_in_bounds = (candidate_h2_index + j) % num_hits;
            const auto h2_index = int(module_pairs[mp_idx-1].hit_start) + index_in_bounds;

            // Check that the phi difference is within tolerance using
            // modulo arithmetic.
            const int16_t phi_diff = phi[h2_index] - extrapolated_phi;
            const int16_t abs_phi_diff = phi_diff < 0 ? -phi_diff : phi_diff;
            if (abs_phi_diff > constants::phi_tolerance_i16)
            {
                break;
            }

            // Extract the hit.
            const auto h2 = hits[h2_index];

            // Work out the scattering.
            const auto dz    = h2.z - h0.z;
            const auto predx = h0.x + tx * dz;
            const auto predy = h0.y + ty * dz;
            const auto dx    = predx - h2.x;
            const auto dy    = predy - h2.y;

            const auto scatter = dx*dx + dy*dy;

            // If this is the best yet, then store it.
            if (scatter < best_fit)
            {
                best_fit = scatter;
                best_h2 = h2_index;
            }
        }

        // Condition for finding a h2.
        if (best_h2 != -1)
        {
            // Mark the hit as used.
            is_used[best_h2] = true;

            // Extend the tracklet to a track.
            if (number_of_hits == 3)
            {
                // Flag the tracklet's hits as used.
                is_used[tracklet->hits[0]] = true;
                is_used[tracklet->hits[1]] = true;
                is_used[tracklet->hits[2]] = true;

                track_number = num_tracks;

                // Store the track.
                tracks[num_tracks].hits[0] = tracklet->hits[0];
                tracks[num_tracks].hits[1] = tracklet->hits[1];
                tracks[num_tracks].hits[2] = tracklet->hits[2];
                tracks[num_tracks].hits[3] = best_h2;
                tracks[num_tracks].num_hits = 4;

                num_tracks++;
            }
            // Add the hit to the existing track.
            else
            {
                track->hits[track->num_hits++] = best_h2;
            }

            // Add to the bag of tracks to follow.
            if (number_of_hits + 1 < constants::max_track_size)
            {
                const auto ttf_p = tracks_to_follow % constants::max_tracks_to_follow;
                track_mask[ttf_p] = track_number;
                tracks_to_follow++;
            }
        }

        // A track just skipped a module, we'll keep it for another round.
        else if (skipped_modules < constants::max_skipped_modules)
        {
            // Create the new mask.
            track_number = ((skipped_modules + 1) << bits::skipped_module_position) |
                            (full_track_number & (bits::seed | bits::track_number));

            // Add to the bag of tracks to follow.
            const auto ttf_p = tracks_to_follow % constants::max_tracks_to_follow;
            track_mask[ttf_p] = track_number;
            tracks_to_follow++;
        }

        // This is a three-hit track.
        else if (number_of_hits == 3)
        {
            three_hit_tracks[num_three_hit_tracks].hits[0] = tracklet->hits[0];
            three_hit_tracks[num_three_hit_tracks].hits[1] = tracklet->hits[1];
            three_hit_tracks[num_three_hit_tracks].hits[2] = tracklet->hits[2];

            num_three_hit_tracks++;
        }
    }
}

std::tuple<int16_t, int16_t> find_forward_candidate(
    ModulePair &module_pair,
    int16_t *phi,
    const Hit &h0,
    const float tx,
    const float ty,
    const float dz)
{
    const auto predx = tx * dz;
    const auto predy = ty * dz;
    const auto x_prediction = h0.x + predx;
    const auto y_prediction = h0.y + predy;
    const auto track_extrapolation_phi = hit_phi_16(x_prediction,
                                                    y_prediction);

    return {binary_search_leftmost(
            phi + int(module_pair.hit_start),
            int(module_pair.num_hits),
            static_cast<int16_t>(track_extrapolation_phi - constants::phi_tolerance_i16)),
            track_extrapolation_phi};
}
