#include "base/base.h"
#include "io/io.h"
#include "io/logging.h"
#include "basis/basis.h"
#include "symmetry/symmetry.h"
#include "integrals/base.h"

#include <chrono>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <format>
#include <iomanip>
#include <sstream>
#include <string>
#include <sstream>

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

namespace fs = std::filesystem;
using SystemClock = std::chrono::system_clock;

static std::string format_time(SystemClock::time_point tp)
{
    const std::time_t t = SystemClock::to_time_t(tp);
    std::tm tm{};

#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif

    std::ostringstream os;
    os << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return os.str();
}

int main(int argc, const char *argv[])
{
    const auto program_start = SystemClock::now();

    logging(LogLevel::Info, "Program Started On :", format_time(program_start));
    logging(LogLevel::Info, "Current Working Directory :", fs::current_path().string());

    if (argc != 2)
    {
        logging(LogLevel::Error, "Usage :", std::format("{} <input file>", argv[0]));
        return EXIT_FAILURE;
    }

    const std::string input_file = argv[1];
    std::ifstream input_stream(input_file);
    if (!input_stream)
    {
        logging(LogLevel::Error, "Input Error :", "Failed to open input file");
        return EXIT_FAILURE;
    }

    // Core objects
    Calculator calculator{};
    Molecule molecule{};

    // Parse input
    if (auto res = read_input(input_stream, calculator, molecule); !res)
    {
        logging(LogLevel::Error, "Input Parsing Failed :", res.error());
        return EXIT_FAILURE;
    }

    logging(LogLevel::Info, "Input Parsing :", "Successful");

    // Detect Symmetry
    if (!calculator.use_pgsymmetry)
    {
        logging(LogLevel::Info, "Symmetry Detection :", "Symmetry detection is turned off by request");
    }

    logging(LogLevel::Info, "Symmetry Detection :", "We use libmsym library to detect point groups");

    if (auto res = detectSymmetry(molecule); !res)
    {
        logging(LogLevel::Error, "Symmetry Detection Failed :", res.error());
        return EXIT_FAILURE;
    }

    logging(LogLevel::Info, "Symmetry Detection :", "Successful");
    logging(LogLevel::Info, "Point Group :", molecule.point_group);

    logging(LogLevel::Info, "Input Coordinates :", "");

    // get input coordinates
    for (std::size_t index = 0; index < molecule.natoms; ++index)
    {
        std::string cstr;
        std::ostringstream astream;
        astream << std::setw(5) << std::right << molecule.atomic_numbers[index];
        cstr += astream.str();

        for (std::size_t cindex = 0; cindex < 3; ++cindex)
        {
            std::ostringstream oss;
            oss << std::setw(10) << std::setprecision(3) << std::fixed << molecule.coordinates[3 * index + cindex];
            cstr += oss.str();
        }
        logging(LogLevel::Info, "", cstr);
    }

    // get reoriented coordinates
    if (molecule.is_reoriented)
    {
        logging(LogLevel::Info, "Standard Coordinates :", "");
        for (std::size_t index = 0; index < molecule.natoms; ++index)
        {
            std::string cstr;
            std::ostringstream astream;
            astream << std::setw(5) << std::right << molecule.atomic_numbers[index];
            cstr += astream.str();

            for (std::size_t cindex = 0; cindex < 3; ++cindex)
            {
                std::ostringstream oss;
                oss << std::setw(10) << std::setprecision(3) << std::fixed << molecule.coordinates_standard[3 * index + cindex];
                cstr += oss.str();
            }
            logging(LogLevel::Info, "", cstr);
        }
    }

    if (calculator.basis_name.empty())
    {
        logging(LogLevel::Error, "Basis Error :", "No basis set file specified");
        return EXIT_FAILURE;
    }

    // Parse basis sets
    const fs::path gbs_path = calculator.basis_path + "/" + calculator.basis_name;
    logging(LogLevel::Info, "Reading Basis Set :", gbs_path.string());

    Basis basis;

    try
    {
        ShellType shell_type = calculator.basis_type.compare("cartesian") == 0 ? ShellType::Cartesian : ShellType::Spherical;
        basis = read_gbs_basis(gbs_path, molecule, shell_type); // cartesian or pure
    }
    catch (const std::exception &e)
    {
        logging(LogLevel::Error, "Basis Parsing Failed :", e.what());
        return EXIT_FAILURE;
    }

    logging(LogLevel::Info, "Basis Construction :", std::format("Generated {} Shells and {} contracted functions", basis.nshells(), basis.nbf()));

    // Find number of shell pairs
    const std::size_t ns = basis.nshells();

    // Store all shell pairs
    std::vector<ShellPair> shell_pairs = {};
    shell_pairs.reserve(ns * (ns + 1) * 0.5);

    for (std::size_t ii = 0; ii < ns; ii++)
    {
        for (std::size_t jj = 0; jj <= ii; jj++)
        {
            shell_pairs.emplace_back(ShellPair{basis.shells[ii], basis.shells[jj]});
        }
    }

    const auto program_end = SystemClock::now();
    const std::chrono::duration<double> elapsed = program_end - program_start;

    logging(LogLevel::Info, "Total Wall Time :", std::format("{:.6f} seconds", elapsed.count()));

    return EXIT_SUCCESS;
}
