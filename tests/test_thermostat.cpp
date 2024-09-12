#include <gtest/gtest.h>
#include "berendsen_thermostat.h"


TEST(BerendsenThermostatTest, VelocityRescaling) {

    // Initialize a simple system with known velocities
    Positions_t positions(3, 2);
    positions.setZero();
    Velocities_t velocities(3, 2);
    velocities << 1.0, 0.0, 0.0, 1.0, 0.0, 0.0;
    Atoms atoms(positions, velocities);

    double target_temperature = 2.0;  // Arbitrary target temperature
    double timestep = 0.001 ;
    double relaxation_time = 0.01;
    int num_steps = 1000;  // Number of steps to apply the thermostat

    // Apply the Berendsen thermostat iteratively
    for (int step = 0; step < num_steps; ++step) {
        berendsen_thermostat(atoms, target_temperature, timestep, relaxation_time);
    }

    // Check if velocities are scaled correctly
    double kinetic_energy_after = 0.5 * atoms.velocities.square().sum();
    double temperature_after = temperature_cur(atoms);

    EXPECT_NEAR(temperature_after, target_temperature, 1e-5);
}
