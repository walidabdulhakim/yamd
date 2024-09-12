//
// Created by walid on 7/3/24.
//

#ifndef YAMD_LJ_DIRECT_SUMMATION_H
#define YAMD_LJ_DIRECT_SUMMATION_H

#include "atoms.h"

double lj_potential(double r, double epsilon, double sigma);
double lj_potential_derivative(double r, double epsilon, double sigma);
double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

#endif // YAMD_LJ_DIRECT_SUMMATION_H
