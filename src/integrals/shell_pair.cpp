#include "shell_pair.h"
#include "math/math.h"

#include <cmath>

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

ShellPair::ShellPair(const Shell &shellA, const Shell &shellB) : shellA(shellA), shellB(shellB), tot_momentumA(shellA.L), tot_momentumB(shellB.L), centerA(shellA.center), centerB(shellB.center)
{
    // Compute distance vector AB
    AB = {
        centerA[0] - centerB[0],
        centerA[1] - centerB[1],
        centerA[2] - centerB[2]};

    // Number of primitive pairs
    const std::size_t na = shellA.exponents.size();
    const std::size_t nb = shellB.exponents.size();

    // Allocate storage for all primitive pairs
    alpha.reserve(na * nb);
    prefac.reserve(na * nb);
    Px.reserve(na * nb);
    Py.reserve(na * nb);
    Pz.reserve(na * nb);

    // Distance squared |AB|**2
    const double AB2 = dot_product(AB, AB);

    // Precompute data for each primitive pair (i,j)
    for (std::size_t i = 0; i < na; ++i)
    {
        const double ai = shellA.exponents[i];
        const double ni = shellA.prim_norms[i];
        const double ci = shellA.coefficients[i];

        for (std::size_t j = 0; j < nb; ++j)
        {
            const double bj = shellB.exponents[j];
            const double nj = shellB.prim_norms[j];
            const double cj = shellB.coefficients[j];

            // 1. Combined exponent: α_ij = α_i + β_j
            const double a = ai + bj;
            alpha.push_back(a);

            // 2. Prefactor includes:
            //    - Gaussian product theorem: exp(-α_i*β_j*|AB|²/(α_i+β_j))
            //    - Contraction coefficients: c_i * d_j
            //    - Primitive normalizations: N_i * N_j
            const double mu = ai * bj / a;
            prefac.push_back(ci * cj * ni * nj * std::exp(-1 * mu * AB2));

            // 3. Gaussian product center P = (α_i * A + β_j * B) / α_ij
            Px.push_back((ai * centerA[0] + bj * centerB[0]) / a);
            Py.push_back((ai * centerA[1] + bj * centerB[1]) / a);
            Pz.push_back((ai * centerA[2] + bj * centerB[2]) / a);
        }
    }
}

std::vector<ShellPair> build_shell_pairs(const Basis &basis)
{
    std::size_t nshells = basis.nshells();
    std::vector<ShellPair> pairs;

    // Reserve space for unique pairs: N*(N+1)/2
    pairs.reserve(nshells * (nshells + 1) / 2);

    // Build only unique pairs (i,j) where i <= j
    for (std::size_t i = 0; i < nshells; ++i)
    {
        for (std::size_t j = i; j < nshells; ++j)
        {
            pairs.emplace_back(basis.shells[i], basis.shells[j]);
        }
    }

    return pairs;
}

std::vector<ShellPair> build_shell_pairs_matrix(const Basis &basis)
{
    std::size_t nshells = basis.nshells();
    std::vector<ShellPair> pairs;

    // Reserve space for all pairs
    pairs.reserve(nshells * nshells);

    for (std::size_t i = 0; i < nshells; ++i)
    {
        for (std::size_t j = 0; j < nshells; ++j)
        {
            pairs.emplace_back(basis.shells[i], basis.shells[j]);
        }
    }

    return pairs;
}