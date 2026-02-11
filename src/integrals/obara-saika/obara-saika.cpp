#include "obara-saika.h"
#include "integrals/shell_pair.h"

#include <cmath>
#include <numbers>

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

static double ObaraSaika::Overlap::computePrimitive1D(int lA, int lB, double PA, double PB, double gamma)
{
    // Base case: S(0,0) = 1
    if (lA == 0 && lB == 0)
    {
        return 1.0;
    }

    // Recursion relations from Obara-Saika scheme
    // S(a+1,b) = PA * S(a,b) + gamma * [a*S(a-1,b) + b*S(a,b-1)]
    // S(a,b+1) = PB * S(a,b) + gamma * [a*S(a-1,b) + b*S(a,b-1)]

    const int max_am = lA + lB;
    std::vector<std::vector<double>> S(max_am + 2, std::vector<double>(max_am + 2, 0.0));

    // Base case
    S[0][0] = 1.0;

    // Build up using recursion
    for (int a = 0; a <= lA; ++a)
    {
        for (int b = 0; b <= lB; ++b)
        {
            if (a == 0 && b == 0)
                continue;

            if (a > 0)
            {
                // Increment a
                S[a][b] = PA * S[a - 1][b];
                if (a > 1)
                    S[a][b] += gamma * a * S[a - 2][b] * 0.5; // Note: factor of 0.5 in gamma definition
                if (b > 0)
                    S[a][b] += gamma * b * S[a - 1][b - 1] * 0.5;
            }
            else
            {
                // Increment b
                S[a][b] = PB * S[a][b - 1];
                if (b > 1)
                    S[a][b] += gamma * b * S[a][b - 2] * 0.5;
            }
        }
    }

    return S[lA][lB];
}

double ObaraSaika::Overlap::computePrimtive3D(const std::array<int, 3> &am_a, const std::array<int, 3> &am_b, const ShellPair &pair, std::size_t prim_idx)
{
    // Extract precomputed data for this primitive pair
    double alpha_ij = pair.alpha[prim_idx]; // α_i + β_j
    double Px = pair.Px[prim_idx];
    double Py = pair.Py[prim_idx];
    double Pz = pair.Pz[prim_idx];

    // γ = 1/(2*α_ij) used in Obara-Saika recursion
    double gamma = 0.5 / alpha_ij;

    // Compute P - A and P - B vectors
    double PAx = Px - pair.centerA[0];
    double PAy = Py - pair.centerA[1];
    double PAz = Pz - pair.centerA[2];

    double PBx = Px - pair.centerB[0];
    double PBy = Py - pair.centerB[1];
    double PBz = Pz - pair.centerB[2];

    // 3D overlap = product of 1D overlaps in x, y, z
    double Sx = ObaraSaika::Overlap::computePrimitive1D(am_a[0], am_b[0], PAx, PBx, gamma);
    double Sy = ObaraSaika::Overlap::computePrimitive1D(am_a[1], am_b[1], PAy, PBy, gamma);
    double Sz = ObaraSaika::Overlap::computePrimitive1D(am_a[2], am_b[2], PAz, PBz, gamma);

    // Include normalization factor from Gaussian product
    // For 3D Gaussians: (pi/α_ij)^(3/2)
    using std::numbers::pi;
    double norm = std::pow(pi / alpha_ij, 1.5);

    return norm * Sx * Sy * Sz;
}

double ObaraSaika::Overlap::computeContracted(const ContractedView &bf_a, const ContractedView &bf_b, const ShellPair &pair)
{
    // Accumulate contribution from all primitive pairs
    double overlap = 0.0;

    std::size_t nprimA = bf_a.exponents().size();
    std::size_t nprimB = bf_b.exponents().size();

    std::size_t prim_idx = 0;
    for (std::size_t i = 0; i < nprimA; ++i)
    {
        for (std::size_t j = 0; j < nprimB; ++j)
        {
            // Compute primitive overlap
            double S_ij = ObaraSaika::Overlap::computePrimtive3D(bf_a.am, bf_b.am, pair, prim_idx);

            // Multiply by contraction coefficients and normalizations (in prefac)
            overlap += pair.prefac[prim_idx] * S_ij;

            ++prim_idx;
        }
    }

    return overlap;
}

std::vector<double> ObaraSaika::Overlap::computeOverlap(const Basis &basis)
{
    std::size_t nbf = basis.nbf();
    std::size_t nshells = basis.nshells();

    // Allocate overlap matrix (nbf × nbf)
    std::vector<double> S(nbf * nbf, 0.0);

    // Build shell pairs (unique pairs only)
    auto shell_pairs = build_shell_pairs(basis);

    // Loop over all contracted basis functions
    std::size_t mu_global = 0; // Global basis function index for shell i
    for (std::size_t ishell = 0; ishell < nshells; ++ishell)
    {
        const auto &shell_i = basis.shells[ishell];

        // Get contracted functions for this shell
        std::size_t nbf_i = 0;
        for (const auto &bf : basis.functions)
        {
            if (bf.shell == &shell_i)
                ++nbf_i;
        }

        std::size_t nu_global = 0; // Global basis function index for shell j
        for (std::size_t jshell = 0; jshell < nshells; ++jshell)
        {
            const auto &shell_j = basis.shells[jshell];

            // Get shell pair
            std::size_t pair_idx = pair_index(ishell, jshell, nshells);
            const auto &pair = shell_pairs[pair_idx];

            // Get contracted functions for shell j
            std::size_t nbf_j = 0;
            for (const auto &bf : basis.functions)
            {
                if (bf.shell == &shell_j)
                    ++nbf_j;
            }

            // Compute overlaps for all function pairs in these shells
            std::size_t mu_local = 0;
            for (const auto &bf_a : basis.functions)
            {
                if (bf_a.shell != &shell_i)
                    continue;

                std::size_t nu_local = 0;
                for (const auto &bf_b : basis.functions)
                {
                    if (bf_b.shell != &shell_j)
                        continue;

                    // Compute S(μ,ν)
                    double overlap = ObaraSaika::Overlap::computeContracted(bf_a, bf_b, pair);

                    // Store in matrix (row-major)
                    std::size_t mu = mu_global + mu_local;
                    std::size_t nu = nu_global + nu_local;
                    S[mu * nbf + nu] = overlap;

                    ++nu_local;
                }
                ++mu_local;
            }

            nu_global += nbf_j;
        }

        mu_global += nbf_i;
    }

    return S;
}