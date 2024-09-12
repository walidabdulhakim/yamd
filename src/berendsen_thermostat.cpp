#include "berendsen_thermostat.h"
#include <cmath>
#include "iostream"

double temperature_cur(const Atoms atoms) {
    return atoms.kinetic_energy() / (KB_EV_K_3_2 * atoms.nb_atoms());
}

void berendsen_thermostat(Atoms &atoms, double target_temperature, double dt, double relaxation_time) {
    double current_temperature = temperature_cur(atoms);

    if (current_temperature > 0.) {
        double ratio = target_temperature / current_temperature;
        double factor = (ratio - 1) * dt / relaxation_time;
        double lambda = sqrt(1 + factor);
        atoms.velocities *= lambda;
    }
}

void berendsen_thermostat_decomposed(Atoms &atoms, double temperature,
                                     double dt, double relaxation_time,
                                     int nb_local, double mass) {
    double temp_cur{atoms.kinetic_energy(nb_local, mass) /
                    (KB_EV_K_3_2 * nb_local)};
    // avoid division by zero, implement berendsen thermostat as described in
    // milestone 5
    double lambda{temp_cur > 0. ? sqrt(1 + (temperature / temp_cur - 1) * dt /
                                               relaxation_time)
                                : 1.};
    atoms.velocities *= lambda;
}

