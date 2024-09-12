//
// Created by walid on 7/3/24.
//

#ifndef YAMD_VERLET_H
#define YAMD_VERLET_H

#include "types.h"

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double timestep);

void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double timestep);

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double dt, Masses_t m);

void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double dt,
                  Masses_t m);

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double dt, double m) ;

void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double dt,
                  double m) ;

#endif // YAMD_VERLET_H
