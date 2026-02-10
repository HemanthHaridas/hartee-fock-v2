#include "base.h"
#include "math/math.h"

#include <cmath>

ShellPair::ShellPair(const Shell &shellA, const Shell &shellB) : shellA(shellA), shellB(shellB), tot_momentumA(shellA.L), tot_momentumB(shellB.L), centerA(shellA.center), centerB(shellB.center)
{
    AB = {centerA[0] - centerB[0], centerA[1] - centerB[1], centerA[2] - centerB[2]};

    const std::size_t na = shellA.exponents.size();
    const std::size_t nb = shellB.exponents.size();

    alpha.reserve(na * nb);
    prefac.reserve(na * nb);
    Px.reserve(na * nb);
    Py.reserve(na * nb);
    Pz.reserve(na * nb);

    const double AB2 = dot_product(AB, AB);

    for (std::size_t i = 0; i < na; ++i)
    {
        const double ai = shellA.exponents[i];
        const double ni = shellA.prim_norms[i];
        const double ci = shellA.coefficients[i];

        for (std::size_t j = 0; j < nb; ++j)
        {
            const double bj = shellB.exponents[j];
            const double nj = shellB.prim_norms[j];
            const double cj = shellB.coefficients[j];

            const double a = ai + bj;
            const double mu = ai * bj / a;

            alpha.push_back(a);

            prefac.push_back(ci * cj * ni * nj * std::exp(-mu * AB2));

            Px.push_back((ai * centerA[0] + bj * centerB[0]) / a);
            Py.push_back((ai * centerA[1] + bj * centerB[1]) / a);
            Pz.push_back((ai * centerA[2] + bj * centerB[2]) / a);
        }
    }
}
