cmake_minimum_required(VERSION 3.15)
project(HartreeFock LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/basis-sets 
    DESTINATION share
)

set(BASIS_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/install/share/basis-sets")

# Collect sources automatically
file(GLOB BASE_SRC
    ${SRC_DIR}/base/*.h
)

file(GLOB BASIS_SRC
    ${SRC_DIR}/basis/*.cpp
    ${SRC_DIR}/basis/*.h
)

file(GLOB IO_SRC
    ${SRC_DIR}/io/*.cpp
    ${SRC_DIR}/io/*.h
)

file(GLOB LOOKUP_SRC
    ${SRC_DIR}/lookup/*.cpp
    ${SRC_DIR}/lookup/*.h
)

file(GLOB SYMM_SRC
    ${SRC_DIR}/symmetry/*.cpp
    ${SRC_DIR}/symmetry/*.h
)

file(GLOB INTERGAL_SRC
    ${SRC_DIR}/integrals/*.cpp
    ${SRC_DIR}/integrals/*.h
)

file(GLOB MAIN_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

configure_file(
    ${SRC_DIR}/base/basis.h.in 
    ${SRC_DIR}/base/basis.h @ONLY
)

# install(
#     FILES ${SRC_DIR}/base/basis.h 
#     DESTINATION include/base
# )

# Create a library from all source files
add_executable(hartree-fock
    ${MAIN_SRC}
    ${BASE_SRC}
    ${BASIS_SRC}
    ${IO_SRC}
    ${LOOKUP_SRC}
    ${INTERGAL_SRC}
    ${SYMM_SRC}
)

# Get the dependencies
set(MSYM_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/external/libmsym/install)
set(MSYM_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/external/libmsym/build)
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install)

include(ExternalProject)

ExternalProject_Add(
    libmsym
    PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/src/external/libmsym
    GIT_REPOSITORY https://github.com/mcodev31/libmsym.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${MSYM_INSTALL_DIR} -DBUILD_SHARED_LIBS:BOOL=OFF -DMSYM_BUILD_EXAMPLES:BOOL=OFF
)

ExternalProject_Get_Property(libmsym install_dir)

# Include directories
# target_include_directories(hartree-fock PUBLIC 
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/base 
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/basis
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/io
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/lookup
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/symmetry
#     ${MSYM_INSTALL_DIR}/include  # Include the external library headers
# )
# Include directories
target_include_directories(hartree-fock PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src # root include path
    ${MSYM_INSTALL_DIR}/include # external lib headers
)

# Link the external library
target_link_libraries(hartree-fock
    ${MSYM_INSTALL_DIR}/lib/libmsym.a
)

# Add dependencies to ensure the external project is built first
add_dependencies(hartree-fock libmsym)

# Install the executable
install(TARGETS hartree-fock DESTINATION)
