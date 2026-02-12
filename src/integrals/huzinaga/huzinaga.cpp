#include "huzinaga.h"
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
* this program.  If not, see <http: //www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/

static double Huzinaga::Overlap::expansionIndex1(int expansionIndex, int lA, int lB, double PA, double PB)
{
    int cMax = std::min(expansionIndex, lA);
    int cMin = std::max(0, expansionIndex - lB);

    double expansionCoeff = 0;
    for (std::size_t ii = cMin; ii <= cMax; ii++)
    {
        auto aux = combination(lA, ii);
        aux = aux * combination(lB, expansionIndex - ii);
        aux = aux * pow(PA, lA - ii);
        aux = aux * pow(PB, lB + ii - expansionIndex);
        expansionCoeff = expansionCoeff + aux;
    }

    return expansionCoeff;
}

static double Huzinaga::Overlap::computePrimitive1D(int lA, int lB, double PA, double PB, double exponentA, double exponentB)
{
    double integral = 0.0;
    double value = 0.0;

    for (std::size_t i = 0; i <= (lA + lB); i++)
    {
        value = double_factorial(2 * i - 1);
        value = value / std::pow(2 * (exponentA + exponentB), i);
        value = value * Huzinaga::Overlap::expansionIndex1(2 * i, lA, lB, PA, PB);
        integral = integral + value;
    }

    return integral;
}

