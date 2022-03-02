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

// Adapted from the LHCb Allen project: https://gitlab.cern.ch/lhcb/Allen

#ifndef _DEFINITIONS_H
#define _DEFINITIONS_H

namespace constants
{
    /// Floating point value of pi.
    static constexpr float pi_f = 3.141592654f;

    /// Maximum scattering.
    static constexpr float max_scatter = 0.08;

    /// Angular tolerance.
    static constexpr float phi_tolerance = 0.045;

    /// Number of module pairs.
    static constexpr unsigned num_module_pairs = 26;

    /// Maximum track size (in hits). (Hit every module pair.)
    static constexpr unsigned max_track_size = num_module_pairs;

    /// Maximum tracks/tracklets to store.
    static constexpr unsigned max_tracks = 280;

    /// Maximum tracks to follow.
    static constexpr unsigned max_tracks_to_follow = 1024;

    /// Maximum number of hits.
    static constexpr unsigned max_hits = max_track_size*max_tracks;

    /// Maximum number of seeding candidates for a track.
    static constexpr unsigned max_seeding_candidates = 8;

    /// Maximum value of 16 bit int as a float.
    static constexpr float max_output_value_i16 = 1 << 15;

    /// Phi conversion factor from float to 16 bit int.
    static constexpr float conversion_factor_i16 = max_output_value_i16 / pi_f;

    /// Phi tolerance as a 16-bit integer.
    static constexpr int16_t phi_tolerance_i16 = static_cast<int16_t>(conversion_factor_i16*phi_tolerance);

    /// Maximum number of skipped modules.
    static constexpr unsigned max_skipped_modules = 1;
}

namespace bits
{
    static constexpr unsigned seed = 0x80000000;
    static constexpr unsigned track_number = 0x0FFFFFFF;
    static constexpr unsigned skipped_modules = 0x70000000;
    static constexpr unsigned skipped_module_position = 28;
}

#endif  /* _DEFINITIONS_H */
