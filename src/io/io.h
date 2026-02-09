#pragma once

#include <istream>
#include <string>
#include <unordered_map>
#include <vector>
#include <expected>

#include "base/base.h"

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

std::expected<SectionMap, std::string> split_into_sections(std::istream &input);
std::expected<Molecule, std::string> parse_geometry(const std::vector<std::string> &lines);
std::expected<Calculator, std::string> parse_calculator(const std::vector<std::string> &lines);
std::expected<void, std::string> read_input(std::istream &input, Calculator &calculator, Molecule &molecule);
