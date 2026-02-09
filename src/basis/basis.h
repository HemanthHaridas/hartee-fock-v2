#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <string>
#include <span>

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

enum class ShellType
{
    Cartesian,
    Spherical
};

// Map shell label ("S", "P", "D", ...) → total angular momentum L
int shell_label_to_L(const std::string &label);

// Generate Cartesian angular momentum tuples for given L
// Example:
//   L = 1 → { {1,0,0}, {0,1,0}, {0,0,1} }
std::vector<std::array<int, 3>> cartesian_shell_order(int L);

// (Optional, later)
// std::vector<std::array<int,3>>
// spherical_shell_order(int L);

// Double factorial: (2n-1)!!
// Required for Gaussian normalization
double double_factorial(int n);

// Primitive normalization constants for a shell of angular momentum L
// Returns vector N_i for each exponent α_i
std::vector<double> primitive_normalization(int L, const std::vector<double> &exponents);
double contraction_normalization(int L, const std::vector<double> &exponents, const std::vector<double> &coefficients, const std::vector<double> &prim_norms);

// Full Cartesian contracted normalization
// (mostly useful for testing / validation)
std::vector<double> cartesian_normalization(const std::array<int, 3> &am, const std::vector<double> &coefficients, const std::vector<double> &exponents);

// Read a Basis Set Exchange .gbs file and build a Basis
Basis read_gbs_basis(const std::string &filename, const Molecule &molecule, ShellType shell_type);
