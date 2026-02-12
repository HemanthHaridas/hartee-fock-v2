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

 namespace Huzinaga
 {
    namespace Overlap
    {
        static double expansionIndex1(int expansionIndex, int lA, int lB, double PA, double PB);
        static double computePrimitive1D(int lA, int lB, double PA, double PB, double exponentA, double exponentB);
    };
 };