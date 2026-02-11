#pragma once

#include <string>
#include <cstdlib>

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

inline std::string get_basis_path()
{
    // Check environment variable override
    const char *env_path = std::getenv("BASIS_PATH");
    if (env_path && *env_path)
    {
        return std::string(env_path);
    }
    // Fallback to compiled-in install path
    return "/Users/hemanthharidas/Desktop/codes/cpp_projects/hartee-fock-v2/install/share/basis-sets";
}
