#pragma once

#include <string>
#include <cstdlib>

inline std::string get_basis_path()
{
    // Check environment variable override
    const char *env_path = std::getenv("BASIS_PATH");
    if (env_path && *env_path)
    {
        return std::string(env_path);
    }
    // Fallback to compiled-in install path
    return "/Users/hemanthharidas/Desktop/codes/cpp_projects/hartee-fock-v2/install/share/basis-sets";
}
