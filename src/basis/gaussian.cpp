#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <array>
#include <vector>
#include <string>
#include <cctype>

#include "base/base.h"
#include "basis.h"
#include "lookup/elements.h"

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

constexpr double ANGSTROM_TO_BOHR = 1.889726124565062;

struct GbsPrimitive
{
    double exponent;
    double coefficient;
};

struct GbsShell
{
    std::string label; // "S", "P", "D", ...
    std::vector<GbsPrimitive> primitives;
};

static inline void normalize_fortran_exponents(std::string &line)
{
    for (char &c : line)
        if (c == 'D' || c == 'd')
            c = 'E';
}

static bool starts_with_alpha(const std::string &line)
{
    for (char c : line)
    {
        if (!std::isspace(static_cast<unsigned char>(c)))
            return std::isalpha(static_cast<unsigned char>(c));
    }
    return false;
}

static bool is_shell_label(const std::string &s)
{
    return s == "S" || s == "P" || s == "D" ||
           s == "F" || s == "G" || s == "H" ||
           s == "SP";
}

using BasisSet = std::unordered_map<std::string, std::vector<GbsShell>>;

static BasisSet read_gbs(std::istream &in)
{
    BasisSet basis;
    std::string line;
    std::string current_element;

    while (std::getline(in, line))
    {
        if (line.empty())
            continue;

        if (line == "****")
        {
            current_element.clear();
            continue;
        }

        if (line.starts_with("!"))
        {
            continue;
        }
        
        /* Element header */
        {
            std::istringstream header(line);
            std::string symbol;
            int charge;

            if ((header >> symbol >> charge) && header.eof())
            {
                element_from_symbol(symbol); // validate
                current_element = symbol;
                basis.try_emplace(symbol);
                continue;
            }
        }

        if (!starts_with_alpha(line))
            throw std::runtime_error("Expected shell header, got: " + line);

        if (current_element.empty())
            throw std::runtime_error("Shell before element header");

        std::istringstream iss(line);
        std::string label;
        std::size_t nprim;
        double scale = 1.0;

        iss >> label >> nprim >> scale;

        if (!iss || !is_shell_label(label))
            throw std::runtime_error("Malformed shell line: " + line);

        if (label == "SP")
        {
            GbsShell s{"S"}, p{"P"};

            for (std::size_t i = 0; i < nprim; ++i)
            {
                std::getline(in, line);
                normalize_fortran_exponents(line);

                std::istringstream prim(line);
                double expn, cs, cp;
                prim >> expn >> cs >> cp;

                s.primitives.push_back({expn, cs * scale});
                p.primitives.push_back({expn, cp * scale});
            }

            basis[current_element].push_back(std::move(s));
            basis[current_element].push_back(std::move(p));
        }
        else
        {
            GbsShell shell{label};

            for (std::size_t i = 0; i < nprim; ++i)
            {
                std::getline(in, line);
                normalize_fortran_exponents(line);

                std::istringstream prim(line);
                double expn, coeff;
                prim >> expn >> coeff;

                shell.primitives.push_back({expn, coeff * scale});
            }

            basis[current_element].push_back(std::move(shell));
        }
    }

    return basis;
}

Basis read_gbs_basis(const std::string &filename, const Molecule &molecule, ShellType shell_type)
{
    std::ifstream file(filename);
    if (!file)
        throw std::runtime_error("Cannot open basis file: " + filename);

    BasisSet gbs = read_gbs(file);
    Basis basis;

    for (std::size_t a = 0; a < molecule.natoms; ++a)
    {
        std::string element = std::string(element_from_z(molecule.atomic_numbers[a]).symbol);

        auto it = gbs.find(element);
        if (it == gbs.end())
            throw std::runtime_error("No basis for element " + element);

        std::array<double, 3> center = {
            molecule.coordinates[3 * a + 0] * ANGSTROM_TO_BOHR,
            molecule.coordinates[3 * a + 1] * ANGSTROM_TO_BOHR,
            molecule.coordinates[3 * a + 2] * ANGSTROM_TO_BOHR};

        for (const GbsShell &gbs_shell : it->second)
        {
            Shell shell;
            shell.center = center;
            shell.L = shell_label_to_L(gbs_shell.label);

            const std::size_t nprim = gbs_shell.primitives.size();
            shell.exponents.reserve(nprim);
            shell.coefficients.reserve(nprim);

            for (const auto &p : gbs_shell.primitives)
            {
                shell.exponents.push_back(p.exponent);
                shell.coefficients.push_back(p.coefficient);
            }

            // perform normalizations
            shell.prim_norms = primitive_normalization(shell.L, shell.exponents);
            const double Nc = contraction_normalization(shell.L, shell.exponents, shell.coefficients, shell.prim_norms);
            for (double &c : shell.coefficients)
                c *= Nc;

            basis.shells.push_back(std::move(shell));
            const Shell *shell_ptr = &basis.shells.back();

            for (auto am : cartesian_shell_order(shell_ptr->L))
            {
                basis.functions.emplace_back(shell_ptr, am);
            }
        }
    }

    return basis;
}
