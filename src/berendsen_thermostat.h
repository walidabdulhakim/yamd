#ifndef BERENDSEN_THERMOSTAT_H
#define BERENDSEN_THERMOSTAT_H

#include "atoms.h"



/// 3/2 times the exact Boltzman constant in units of eV/K to 9 decimal places.
const double KB_EV_K_3_2{1.5 * KB_EV_K};

double temperature_cur( const Atoms atoms);
void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep, double relaxation_time);
void berendsen_thermostat_decomposed(Atoms &atoms, double temperature,
                                     double dt, double relaxation_time,
                                     int nb_local, double mass);


#endif // BERENDSEN_THERMOSTAT_H
