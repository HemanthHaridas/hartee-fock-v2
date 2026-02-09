#pragma once

#include <array>
#include <string_view>
#include <cstdint>

/**
 * @brief Element metadata container
 *
 * This is the single authoritative source for periodic table data
 * used throughout the codebase.
 */
struct ElementData
{
    std::string_view symbol;   ///< Chemical symbol (e.g. "C")
    std::uint64_t    Z;        ///< Atomic number
    double           mass;     ///< Atomic mass (amu)
    double           radius;   ///< Covalent radius (Å)
};

/**
 * @brief Periodic table data (H → Es)
 *
 * Data source:
 *   - Atomic masses: PubChem
 *   - Radii: covalent radii (used for heuristics only)
 *
 * This table is immutable and defined in elements.cpp.
 */
extern const std::array<ElementData, 99> periodic_table;

/**
 * @brief Lookup element data by chemical symbol
 *
 * @param symbol Element symbol (e.g. "C")
 * @return Reference to ElementData
 * @throws std::runtime_error if symbol is unknown
 */
const ElementData& element_from_symbol(std::string_view symbol);

/**
 * @brief Lookup element data by atomic number
 *
 * @param Z Atomic number
 * @return Reference to ElementData
 * @throws std::runtime_error if Z is invalid
 */
const ElementData& element_from_z(std::uint64_t Z);
