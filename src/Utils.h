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

#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>

#include "Definitions.h"

inline float signselect(const float& s, const float& a, const float& b)
{
    return (s > 0) ? a : b;
}

float fast_atan2f(const float &y, const float &x)
{
    const float c1 = .25f * constants::pi_f;
    const float c2 = .75f * constants::pi_f;
    const float abs_y = fabsf(y);

    const float x_plus_y = x + abs_y;
    const float x_sub_y = x - abs_y;
    const float y_sub_x = abs_y - x;

    const float nom = signselect(x, x_sub_y, x_plus_y);
    const float den = signselect(x, x_plus_y, y_sub_x);

    const float r = nom / den;
    const float angle = signselect(x, c1, c2) - c1 * r;

    return copysignf(angle, y);
}

int16_t hit_phi_16(const float x, const float y)
{
    const float float_value = fast_atan2f(y, x) * constants::conversion_factor_i16;
    const int16_t int16_value = static_cast<int16_t>(float_value);
    return int16_value;
}

int binary_search_leftmost(const int16_t *array, const unsigned array_size, const int16_t &value)
{
    int l = 0;
    int r = array_size;
    while (l < r)
    {
        const int m = (l + r) / 2;
        const auto array_element = array[m];
        if (value > array_element)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    return l;
}

#endif  /* _UTILS_H */
