#pragma once

#include <cstddef>   
#include <algorithm> 

#include "base/base.h"

/*-----------------------------------------------------------------------------
 * Planck
 * Copyright (C) 2024 Hemanth Haridas, University of Utah
 * Contact: hemanthhari23@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or a later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 ----------------------------------------------------------------------------*/

struct ShellPair
{
    // pointers to Shells
    const Shell &shellA;
    const Shell &shellB;

    // L values
    int tot_momentumA;
    int tot_momentumB;

    // Centers
    std::array<double, 3> centerA;
    std::array<double, 3> centerB;

    // Distance vector AB
    std::array<double, 3> AB;

    // Precomputed primitive-pair data
    std::vector<double> alpha; // α_i + β_j
    std::vector<double> prefac;
    std::vector<double> Px, Py, Pz;

    ShellPair(const Shell &shellA, const Shell &shellB);
};

std::vector<ShellPair> build_shell_pairs(const Basis &basis);
std::vector<ShellPair> build_shell_pairs_matrix(const Basis &basis);

// Compute the shell pair index for (i,j) given total number of shells
inline std::size_t pair_index(std::size_t i, std::size_t j, std::size_t nshells)
{
    if (i > j)
        std::swap(i, j);

    // Formula for mapping (i,j) to unique pair index
    return i * nshells - i * (i - 1) / 2 + (j - i);
}
