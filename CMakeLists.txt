cmake_minimum_required(VERSION 3.15)
project(HartreeFock LANGUAGES CXX)

# C++ standard settings
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Source directory
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# Install basis sets
install(
    DIRECTORY ${CMAKE_SOURCE_DIR}/basis-sets 
    DESTINATION share
)

set(BASIS_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/share/basis-sets")

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
    ${SRC_DIR}/integrals/*/*.cpp
    ${SRC_DIR}/integrals/*/*.h
)

file(GLOB MAIN_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

# Configure generated header
configure_file(
    ${SRC_DIR}/base/basis.h.in 
    ${SRC_DIR}/base/basis.h @ONLY
)

# Executable target
add_executable(hartree-fock
    ${MAIN_SRC}
    ${BASE_SRC}
    ${BASIS_SRC}
    ${IO_SRC}
    ${LOOKUP_SRC}
    ${INTEGRAL_SRC}
    ${SYMM_SRC}
)

# Optimization flags per platform
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(STATUS "Configuring for Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ftree-vectorize")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    message(STATUS "Configuring for macOS")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math")
elseif (CMAKE_SYSTEM_NAME STREQUAL "Windows")
    message(STATUS "Configuring for Windows")
    # MSVC uses different flags
    if (MSVC)
        add_compile_options(/O2 /arch:AVX2)
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ftree-vectorize")
    endif()
endif()

# External project: libmsym
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
target_include_directories(hartree-fock PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${MSYM_INSTALL_DIR}/include
)

# Link external library
target_link_libraries(hartree-fock
    ${MSYM_INSTALL_DIR}/lib/libmsym.a
)

# Ensure libmsym builds before hartree-fock
add_dependencies(hartree-fock libmsym)

# Install executable
install(TARGETS hartree-fock DESTINATION bin)
