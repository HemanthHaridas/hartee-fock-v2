#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>
#include <span>

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

struct Molecule
{
    std::size_t natoms = 0;

    std::vector<std::uint64_t> atomic_numbers; // Z
    std::vector<double> atomic_masses;         // amu

    // Cartesian coordinates in BOHR
    // Layout: [x0, y0, z0, x1, y1, z1, ...]
    std::vector<double> coordinates;
    std::vector<double> coordinates_standard;

    std::string point_group;
    bool is_reoriented = false;

    void clear() noexcept
    {
        natoms = 0;
        atomic_numbers.clear();
        atomic_masses.clear();
        coordinates.clear();
    }
};

struct Shell
{
    // Center in BOHR
    std::array<double, 3> center{};

    // Total angular momentum (L = lx + ly + lz)
    int L = 0;

    // Primitive data (SoA layout)
    std::vector<double> exponents;    // α_i
    std::vector<double> coefficients; // c_i
    std::vector<double> prim_norms;   // primitive normalization constants

    std::size_t nprimitives() const noexcept
    {
        return exponents.size();
    }
};

struct ContractedView
{
    // Pointer to owning shell (non-owning, stable lifetime)
    const Shell *shell = nullptr;

    // Cartesian angular momentum component (lx, ly, lz)
    std::array<int, 3> am{};

    // Convenience accessors
    std::span<const double> exponents() const noexcept
    {
        return shell->exponents;
    }

    std::span<const double> coefficients() const noexcept
    {
        return shell->coefficients;
    }

    std::span<const double> primitive_norms() const noexcept
    {
        return shell->prim_norms;
    }

    const std::array<double, 3> &center() const noexcept
    {
        return shell->center;
    }
};

struct Basis
{
    // Owns shells
    std::vector<Shell> shells;

    // Lightweight contracted-function views
    std::vector<ContractedView> functions;

    std::size_t nshells() const noexcept
    {
        return shells.size();
    }

    std::size_t nbf() const noexcept
    {
        return functions.size();
    }

    void clear()
    {
        shells.clear();
        functions.clear();
    }
};

enum IntegralEngine
{
    MD,
    THO,
    OS
};

struct Calculator
{
    // Input
    std::string basis_name;
    std::string basis_path;
    std::string method;
    std::string calc_type;  // calculation type
    std::string coord_type; // cartesian / z-matrix
    std::string basis_type; // cartesian / spherical

    IntegralEngine integral_engine = IntegralEngine::OS;

    int max_iter = 50;
    int max_scf = 50;
    int diis_dim = 10;

    double tol_scf = 1e-10;
    double tol_eri = 1e-10;

    bool use_pgsymmetry = true;
    bool use_diis = true;

    int charge = 0;
    int multiplicity = 1;
    int tot_electrons = 0;
    
    // Output
    double final_energy = 0.0;
    bool converged = false;

    // SCF matrices (row-major)
    std::vector<double> C; // MO coefficients
    std::vector<double> D; // density matrix

    void resize(std::size_t nbf)
    {
        C.assign(nbf * nbf, 0.0);
        D.assign(nbf * nbf, 0.0);
    }

    void reset()
    {
        final_energy = 0.0;
        converged = false;
        C.clear();
        D.clear();
    }
};
