#pragma once

#include "base/base.h"

struct ShellPair
{
    // pointers to Shells
    const Shell &shellA;
    const Shell &shellB;

    // L values
    int tot_momentumA;
    int tot_momentumB;

    // Centers
    std::array<double, 3> centerA;
    std::array<double, 3> centerB;

    // Distance vector AB
    std::array<double, 3> AB;

    // Precomputed primitive-pair data
    std::vector<double> alpha; // α_i + β_j
    std::vector<double> prefac;
    std::vector<double> Px, Py, Pz;

    ShellPair(const Shell &shellA, const Shell &shellB);
};