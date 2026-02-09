#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <expected>
#include <sstream>
#include <cctype>

#include "io.h"
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

using SectionMap = std::unordered_map<std::string, std::vector<std::string>>;

static void trim(std::string &s)
{
    const auto first = s.find_first_not_of(" \t");
    if (first == std::string::npos)
    {
        s.clear();
        return;
    }

    const auto last = s.find_last_not_of(" \t");
    s = s.substr(first, last - first + 1);
}

std::expected<SectionMap, std::string>split_into_sections(std::istream &input)
{
    SectionMap sections;

    std::string line;
    std::string current;
    bool in_section = false;

    while (std::getline(input, line))
    {
        trim(line);

        if (line.empty() || line.front() == '#')
            continue;

        // Section header?
        if (line.front() == '[' && line.back() == ']')
        {
            const std::string tag = line.substr(1, line.size() - 2);

            // END tag
            if (tag.starts_with("END "))
            {
                const std::string end_name = tag.substr(4);

                if (!in_section)
                    return std::unexpected("END without active section: " + end_name);

                if (end_name != current)
                    return std::unexpected(
                        "Mismatched END section. Expected END " + current +
                        ", got END " + end_name);

                current.clear();
                in_section = false;
                continue;
            }

            // START tag
            if (in_section)
                return std::unexpected(
                    "Nested section [" + tag + "] inside [" + current + "]");

            current = tag;
            in_section = true;
            sections[current]; // create entry
            continue;
        }

        if (in_section)
        {
            sections[current].push_back(line);
        }
    }

    if (in_section)
        return std::unexpected("Unterminated section: " + current);

    if (sections.empty())
        return std::unexpected("No sections found in input");

    return sections;
}

std::expected<Molecule, std::string>parse_geometry(const std::vector<std::string> &lines)
{
    if (lines.empty())
        return std::unexpected("Empty GEOM section");

    std::size_t natoms = 0;
    {
        std::istringstream iss(lines[0]);
        if (!(iss >> natoms) || natoms == 0)
            return std::unexpected("Invalid atom count in GEOM section");
    }

    if (lines.size() != natoms + 1)
        return std::unexpected("GEOM atom count does not match number of lines");

    Molecule mol;
    mol.natoms = natoms;
    mol.atomic_numbers.resize(natoms);
    mol.atomic_masses.resize(natoms);
    mol.coordinates.resize(3 * natoms);

    for (std::size_t i = 0; i < natoms; ++i)
    {
        std::istringstream iss(lines[i + 1]);

        std::string symbol;
        double x, y, z;

        if (!(iss >> symbol >> x >> y >> z))
            return std::unexpected("Malformed GEOM line: " + lines[i + 1]);

        try
        {
            const auto &el = element_from_symbol(symbol);
            mol.atomic_numbers[i] = el.Z;
            mol.atomic_masses[i] = el.mass;
        }
        catch (...)
        {
            return std::unexpected("Unknown atomic symbol: " + symbol);
        }

        mol.coordinates[3 * i + 0] = x;
        mol.coordinates[3 * i + 1] = y;
        mol.coordinates[3 * i + 2] = z;
    }

    return mol;
}

std::expected<Calculator, std::string>parse_calculator(const std::vector<std::string> &lines)
{
    if (lines.size() < 2)
        return std::unexpected("Incomplete CALC section");

    Calculator calc;

    {
        std::istringstream iss(lines[0]);
        if (!(iss >> calc.basis_name >> calc.method >> calc.basis_type))
            return std::unexpected("Malformed CALC header: " + lines[0]);
    }

    {
        std::istringstream iss(lines[1]);
        if (!(iss >> calc.charge >> calc.multiplicity))
            return std::unexpected("Malformed charge/multiplicity line: " + lines[1]);
    }

    if (calc.multiplicity <= 0)
        return std::unexpected("Invalid spin multiplicity");

    return calc;
}

std::expected<void, std::string>read_input(std::istream &input, Calculator &calc, Molecule &mol)
{
    if (!input)
        return std::unexpected("Invalid input stream");

    auto sections = split_into_sections(input);
    if (!sections)
        return std::unexpected(sections.error());

    // GEOM
    auto geom_it = sections->find("GEOM");
    if (geom_it == sections->end())
        return std::unexpected("Missing required [GEOM] section");

    auto geom = parse_geometry(geom_it->second);
    if (!geom)
        return std::unexpected(geom.error());

    // CALC
    auto calc_it = sections->find("CALC");
    if (calc_it == sections->end())
        return std::unexpected("Missing required [CALC] section");

    auto calc_parsed = parse_calculator(calc_it->second);
    if (!calc_parsed)
        return std::unexpected(calc_parsed.error());

    mol  = std::move(*geom);
    calc = std::move(*calc_parsed);

    return {};
}
