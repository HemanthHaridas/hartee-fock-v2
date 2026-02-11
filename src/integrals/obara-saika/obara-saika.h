#pragma once

#include "base/base.h"
#include "integrals/shell_pair.h"

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

namespace ObaraSaika
{
    namespace Overlap
    {
        static double computePrimitive1D(int lA, int lB, double PA, double PB, double gamma);
        double computePrimtive3D(const std::array<int, 3> &am_a, const std::array<int, 3> &am_b, const ShellPair &pair, std::size_t prim_idx);
        double computeContracted(const ContractedView &bf_a, const ContractedView &bf_b, const ShellPair &pair);
        std::vector<double> computeOverlap(const Basis &basis);
    };

    namespace Kinetic
    {

    };
};