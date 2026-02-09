#include <algorithm>
#include <iostream>
#include <iomanip>
#include <set>
#include <string.h>

#include "symmetry.h"
#include "wrapper.h"

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

// std::expected<void, std::string> detectSymmetry(Molecule &molecule)
// {
//     msym_thresholds_t sloppy_thresholds = {
//         0.08,   // zero
//         0.1,    // geometry
//         0.1,    // angle
//         0.06,   // equivalence
//         1.0e-1, // permutation
//         1.0e-3, // eigfact
//         0.1     // orthogonalization
//     };

//     // Need to create a msym_element_arry for using libmsym
//     // this section of the code is unfortunately in C
//     msym_error_t ret = MSYM_SUCCESS;
//     msym_element_t *elements = NULL;
//     char point_group[10];
//     unsigned int length = molecule.natoms;

//     msym_element_t *geometry = NULL;
//     msym_element_t *atom;
//     atom = (msym_element_t *)malloc(molecule.natoms * sizeof(msym_element_t));
//     memset(atom, 0, molecule.natoms * sizeof(msym_element_t));

//     for (std::uint64_t ii = 0; ii < molecule.natoms; ii++)
//     {
//         atom[ii].m = molecule.atomic_masses[ii];  // mass
//         atom[ii].n = molecule.atomic_numbers[ii]; // nuclear charge

//         // now set positions
//         atom[ii].v[0] = molecule.coordinates[ii * 3 + 0];
//         atom[ii].v[1] = molecule.coordinates[ii * 3 + 1];
//         atom[ii].v[2] = molecule.coordinates[ii * 3 + 2];
//     }
//     geometry = atom;
//     msym_context symm_context = msymCreateContext();
//     msymSetThresholds(symm_context, &sloppy_thresholds);

//     // need to initialize the eigen matrix to hold standardized coordinates
//     molecule.coordinates_standard.resize(molecule.natoms * 3);

//     // now check the point group of the molecule
//     if (MSYM_SUCCESS != (ret = msymSetElements(symm_context, length, geometry)))
//     {
//         free(geometry);
//         return std::unexpected("Unable to set elements.");
//     }

//     if (MSYM_SUCCESS != (ret = msymFindSymmetry(symm_context)))
//     {
//         free(geometry);
//         molecule.point_group = "C1";
//         molecule.coordinates_standard = molecule.coordinates;
//         molecule.is_reoriented = false;
//     }

//     if (MSYM_SUCCESS != (ret = msymGetPointGroupName(symm_context, sizeof(point_group), point_group)))
//     {
//         free(geometry);
//     }

//     molecule.point_group = point_group;
//     free(geometry);

//     // if it has a Cinf axis, replace 0 with inf
//     if (point_group[1] == '0')
//     {
//         molecule.point_group = molecule.point_group.replace(1, 1, "inf");
//     }

//     // symmetrize the molecule before that
//     msym_element_t *new_geometry = NULL;
//     double symm_error = 0.0;
//     int new_n_atoms = 0;

//     // first check if the symmetrization is success
//     if (MSYM_SUCCESS != (ret = msymSymmetrizeElements(symm_context, &symm_error)))
//     {
//         return std::unexpected("Unable to symmetrize the molecule to a higher point group.");
//     }

//     if (MSYM_SUCCESS != (ret = msymGetElements(symm_context, &new_n_atoms, &new_geometry)))
//     {
//         return std::unexpected("Unable to get the symmetry elements for " + molecule.point_group + "after symmetrization.");
//     }

//     // align the axes correctly
//     if (MSYM_SUCCESS != (ret = msymAlignAxes(symm_context)))
//     {
//         return std::unexpected("Unable to align the symmetry axes.");
//     }

//     // now update the new geometry
//     for (std::uint64_t ii = 0; ii < molecule.natoms; ii++)
//     {
//         molecule.coordinates_standard[ii * 3 + 0] = (new_geometry[ii].v[0]);
//         molecule.coordinates_standard[ii * 3 + 1] = (new_geometry[ii].v[1]);
//         molecule.coordinates_standard[ii * 3 + 2] = (new_geometry[ii].v[2]);
//     }
//     molecule.is_reoriented = true;
//     msymReleaseContext(symm_context);
//     // free(atom);
//     return {};
// }

// Full implementation of detectSymmetry
std::expected<void, std::string> detectSymmetry(Molecule &molecule)
{
    try
    {
        SymmetryContext ctx;

        msym_thresholds_t sloppy_thresholds = {
            0.08,   // zero
            0.1,    // geometry
            0.1,    // angle
            0.06,   // equivalence
            1.0e-1, // permutation
            1.0e-3, // eigfact
            0.1     // orthogonalization
        };

        msymSetThresholds(ctx.get(), &sloppy_thresholds);

        SymmetryElements atoms(molecule.natoms);
        for (size_t i = 0; i < molecule.natoms; ++i)
        {
            atoms.data()[i].m = molecule.atomic_masses[i];
            atoms.data()[i].n = molecule.atomic_numbers[i];
            atoms.data()[i].v[0] = molecule.coordinates[i * 3 + 0];
            atoms.data()[i].v[1] = molecule.coordinates[i * 3 + 1];
            atoms.data()[i].v[2] = molecule.coordinates[i * 3 + 2];
        }

        if (MSYM_SUCCESS != msymSetElements(ctx.get(), atoms.size(), atoms.data()))
        {
            return std::unexpected("Unable to set elements.");
        }

        if (MSYM_SUCCESS != msymFindSymmetry(ctx.get()))
        {
            molecule.point_group = "C1";
            molecule.coordinates_standard = molecule.coordinates;
            molecule.is_reoriented = false;
            return {};
        }

        char point_group[32];
        if (MSYM_SUCCESS != msymGetPointGroupName(ctx.get(), sizeof(point_group), point_group))
        {
            return std::unexpected("Unable to get point group name.");
        }
        molecule.point_group = point_group;

        if (point_group[1] == '0')
        {
            molecule.point_group.replace(1, 1, "inf");
        }

        double symm_error = 0.0;
        if (MSYM_SUCCESS != msymSymmetrizeElements(ctx.get(), &symm_error))
        {
            return std::unexpected("Unable to symmetrize the molecule.");
        }

        int new_n_atoms = 0;
        msym_element_t *new_geometry = nullptr;
        if (MSYM_SUCCESS != msymGetElements(ctx.get(), &new_n_atoms, &new_geometry))
        {
            return std::unexpected("Unable to get symmetry elements.");
        }

        if (MSYM_SUCCESS != msymAlignAxes(ctx.get()))
        {
            return std::unexpected("Unable to align symmetry axes.");
        }

        molecule.coordinates_standard.resize(molecule.natoms * 3);
        for (size_t i = 0; i < molecule.natoms; ++i)
        {
            molecule.coordinates_standard[i * 3 + 0] = new_geometry[i].v[0];
            molecule.coordinates_standard[i * 3 + 1] = new_geometry[i].v[1];
            molecule.coordinates_standard[i * 3 + 2] = new_geometry[i].v[2];
        }
        molecule.is_reoriented = true;

        return {};
    }
    catch (const std::exception &e)
    {
        return std::unexpected(e.what());
    }
}
