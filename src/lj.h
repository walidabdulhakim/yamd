//
// Created by walid on 7/3/24.
//

#ifndef YAMD_LJ
#define YAMD_LJ

#include "atoms.h"
#include "neighbors.h"

double lj_cutoff_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0, double cutoff=5.0);

#endif // YAMD_LJ
