#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <functional>
#include <expected>
#include <algorithm>
#include <sstream>
#include <cctype>

#include "base/basis.h"
#include "io.h"
#include "lookup/elements.h"
#include <numeric>

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

std::string toLower(const std::string &parsedString)
{
    std::string lowerString = parsedString; // Create a copy to preserve the original string
    std::transform(lowerString.begin(), lowerString.end(), lowerString.begin(), ::tolower);
    return lowerString; // Return the transformed string
}

bool stringToBool(const std::string &parsedString)
{
    std::string upperStr = parsedString;
    std::transform(upperStr.begin(), upperStr.end(), upperStr.begin(), ::toupper); // Convert to uppercase

    if (upperStr == "ON")
    {
        return true; // "ON" maps to true
    }
    else if (upperStr == "OFF")
    {
        return false; // "OFF" maps to false
    }
    else
    {
        throw std::invalid_argument("Invalid string for boolean conversion."); // Handle unexpected input
    }
}

std::expected<void, std::string> check_charge_multiplicity(const Molecule &molecule, const Calculator &calc)
{
    int n_elec = std::accumulate(molecule.atomic_numbers.begin(), molecule.atomic_numbers.end(), 0) - calc.charge;

    if (n_elec <= 0)
        return std::unexpected("Invalid electron count: " + std::to_string(n_elec));

    if (calc.multiplicity <= 0)
        return std::unexpected("Invalid spin multiplicity: " + std::to_string(calc.multiplicity));

    // Parity check: odd multiplicity → even electrons, even multiplicity → odd electrons
    bool parity_ok = ((calc.multiplicity % 2 == 1 && n_elec % 2 == 0) || (calc.multiplicity % 2 == 0 && n_elec % 2 == 1));

    if (!parity_ok)
    {
        return std::unexpected(
            "Parity mismatch: electron count (" + std::to_string(n_elec) + ") incompatible with multiplicity (" + std::to_string(calc.multiplicity) + ")");
    }

    // Max multiplicity check: multiplicity cannot exceed n_elec + 1
    if (calc.multiplicity > n_elec + 1)
    {
        return std::unexpected("Multiplicity (" + std::to_string(calc.multiplicity) + ") exceeds maximum possible for " + std::to_string(n_elec) + " electrons");
    }

    // check if theory and multiplcity match
    if (calc.multiplicity > 1 && calc.method.compare("uhf"))
    {
        return std::unexpected("Multiplicity (" + std::to_string(calc.multiplicity) + ") is incompatibile with " + calc.method);
    }

    return {}; // success
}

std::expected<SectionMap, std::string> split_into_sections(std::istream &input)
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

std::expected<Molecule, std::string> parse_geometry(const std::vector<std::string> &lines)
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

std::expected<Calculator, std::string> parse_calculator(const std::vector<std::string> &lines)
{
    if (lines.size() < 2)
    {
        return std::unexpected("Incomplete CALC section");
    }

    Calculator calc;

    // maps to varaibles
    std::unordered_map<std::string, std::function<void(std::string)>> handlers_setup = {
        // calculation information
        {"CALC_TYPE",   [&calc](std::string value){ calc.calc_type  = toLower(value); }},
        {"THEORY",      [&calc](std::string value){ calc.method     = toLower(value); }},
        {"BASIS",       [&calc](std::string value){ calc.basis_name = toLower(value); }},
        // {"BASIS_PATH",  [&calc](std::string value){ calc.basis_path = toLower(value); }},

        // diis and symmetry information
        {"USE_SYMM",    [&calc](std::string value){ calc.use_pgsymmetry = stringToBool(value); }},
        {"USE_DIIS",    [&calc](std::string value){ calc.use_diis       = stringToBool(value); }},

        // max cycles, charge and multiplicity
        {"MAXITER",     [&calc](std::string value){ calc.max_iter       = std::stoi(value); }},
        {"MAXSCF",      [&calc](std::string value){ calc.max_scf        = std::stoi(value); }},
        {"CHARGE",      [&calc](std::string value){ calc.charge         = std::stoi(value); }},
        {"MULTI",       [&calc](std::string value){ calc.multiplicity   = std::stoi(value); }},
        {"DIIS_DIM",    [&calc](std::string value){ calc.diis_dim       = std::stoi(value); }},

        // tolerances
        {"TOLSCF",      [&calc](std::string value){ calc.tol_scf   = std::stod(value); }},
        {"TOLERI",      [&calc](std::string value){ calc.tol_eri   = std::stod(value); }}
        };

    for (auto line : lines)
    {
        std::istringstream iss(line);
        std::string key, value;

        if (!(iss >> key >> value))
        {
            return std::unexpected("Malformed Input line");
        }

        // search for the key in map
        auto it = handlers_setup.find(key);

        if (it == handlers_setup.end())
        {
            return std::unexpected("Key not found. Check [CALC] block");
        }

        // now set the value
        it->second(value);
    }

    if (calc.multiplicity <= 0)
        return std::unexpected("Invalid spin multiplicity");

    calc.basis_path = get_basis_path();
    return calc;
}

std::expected<void, std::string> read_input(std::istream &input, Calculator &calc, Molecule &mol)
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

    // check if charge and multiplicity match
    if (auto checks = check_charge_multiplicity(*geom, *calc_parsed); !checks)
    {
        return std::unexpected(checks.error());
    }

    mol = std::move(*geom);
    calc = std::move(*calc_parsed);

    return {};
}
