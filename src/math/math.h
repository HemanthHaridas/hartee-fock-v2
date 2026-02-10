#pragma once

#include <vector>
#include <array>
#include <stdexcept>

// Generic dot product for vectors of any numeric type
template <typename T>
T dot_product(const std::vector<T> &a, const std::vector<T> &b)
{
    if (a.size() != b.size())
    {
        throw std::invalid_argument("Vectors must be the same length");
    }
    T result = T{};
    for (size_t i = 0; i < a.size(); ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

// Generic dot product for arrays of any dimension
template <typename T, std::size_t N>
T dot_product(const std::array<T, N>& a, const std::array<T, N>& b)
{
    T result = T{};
    for (std::size_t i = 0; i < N; ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}
