#include "verlet.h"

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double timestep) {

    positions += velocities * timestep + 0.5 * forces * timestep * timestep;
    velocities += 0.5 * forces * timestep;
}

void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double timestep) {
    velocities += 0.5 * forces * timestep;
}



void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double dt, Masses_t m) {
    velocities += (forces.rowwise() / m.transpose()) * dt * .5;
    positions += velocities * dt;
}


void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double dt,
                  Masses_t m) {
    velocities += (forces.rowwise() / m.transpose()) * dt * .5;
}


void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, double dt, double m) {
    velocities += forces / m * dt * .5;
    positions += velocities * dt;
}

void verlet_step2(Velocities_t &velocities, const Forces_t &forces, double dt,
                  double m) {
    velocities += forces / m * dt * .5;
}