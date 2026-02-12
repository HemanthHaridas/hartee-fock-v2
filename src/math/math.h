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
T dot_product(const std::array<T, N> &a, const std::array<T, N> &b)
{
    T result = T{};
    for (std::size_t i = 0; i < N; ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

double double_factorial(int n)
{
    // Base case when n < -1
    if (n < -1)
    {
        return 0;
    }

    // 0!! = 1!! = (-1)!! = 1 by definition
    if (n <= 1)
    {
        return 1;
    }

    double result = 1;
    for (size_t i = n; i > 1; i -= 2)
    {
        result = result * i;
    }

    return result;
}

long double combination(int n, int r)
{
    if (r < 0 || r > n)
        throw std::runtime_error("Invalid r for combination");
    if (r == 0 || r == n)
        return 1.0L;

    // Use symmetry to reduce iterations
    if (r > n - r)
        r = n - r;
    long double result = 1.0L;
    for (int k = 1; k <= r; ++k)
    {
        // Multiply first, then divide — keeps result exact
        result *= (n - r + k);
        result /= k;
    }
    return result;
}