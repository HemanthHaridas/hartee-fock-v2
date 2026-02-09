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

#include "basis.h"

int shell_label_to_L(const std::string &label)
{
    if (label == "S")
        return 0;
    if (label == "P")
        return 1;
    if (label == "D")
        return 2;
    if (label == "F")
        return 3;
    if (label == "G")
        return 4;
    if (label == "H")
        return 5;

    throw std::runtime_error("Unknown shell label: " + label);
}

std::vector<std::array<int, 3>> cartesian_shell_order(int L)
{
    std::vector<std::array<int, 3>> result;
    result.reserve((L + 1) * (L + 2) / 2);

    for (int lx = L; lx >= 0; --lx)
    {
        for (int ly = L - lx; ly >= 0; --ly)
        {
            int lz = L - lx - ly;
            result.push_back({lx, ly, lz});
        }
    }

    return result;
}

std::vector<double> primitive_normalization(int L, const std::vector<double> &exponents)
{
    constexpr double pi = 3.1415926535897932384626433832795;

    const double prefactor = std::pow(2.0, 2.0 * L + 1.5) / std::pow(pi, 1.5);

    std::vector<double> norms;
    norms.reserve(exponents.size());

    for (double alpha : exponents)
    {
        const double n = std::sqrt(std::pow(alpha, L + 1.5) * prefactor);
        norms.push_back(n);
    }

    return norms;
}

double contraction_normalization(int L, const std::vector<double> &exponents, const std::vector<double> &coefficients, const std::vector<double> &prim_norms)
{
    constexpr double pi = 3.1415926535897932384626433832795;

    const std::size_t n = exponents.size();

    if (coefficients.size() != n || prim_norms.size() != n)
        throw std::runtime_error("contraction_normalization: size mismatch");

    double sum = 0.0;

    for (std::size_t i = 0; i < n; ++i)
    {
        const double ai = exponents[i];
        const double ci = coefficients[i];
        const double Ni = prim_norms[i];

        for (std::size_t j = 0; j < n; ++j)
        {
            const double aj = exponents[j];
            const double cj = coefficients[j];
            const double Nj = prim_norms[j];

            const double aij = ai + aj;

            const double Sij =
                std::pow(pi / aij, 1.5) *
                std::pow(2.0 * std::sqrt(ai * aj) / aij, L);

            sum += ci * cj * Ni * Nj * Sij;
        }
    }

    return 1.0 / std::sqrt(sum);
}
