#include "gtest/gtest.h"
#include "verlet.h"

TEST(VerletTest, MultipleAtoms) {
    double dt = 0.01;
    int nb_atoms = 3;
    int steps = 100;

    Positions_t positions(3, nb_atoms);
    positions.setZero();

    Velocities_t velocities(3, nb_atoms);
    velocities.setZero();

    Forces_t forces(3, nb_atoms);
    forces.setOnes();

    for (int i = 0; i < steps; ++i) {
        verlet_step1(positions, velocities, forces, dt);
        verlet_step2(velocities, forces, dt);
    }

    // Correctly calculate expected positions and velocities
    double time = steps * dt;
    Eigen::Vector3d expected_position = 0.5 * forces.col(0) * time * time;
    Eigen::Vector3d expected_velocity = forces.col(0) * time;

    for (int i = 0; i < nb_atoms; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(positions(j, i), expected_position(j), 1e-2);
            EXPECT_NEAR(velocities(j, i), expected_velocity(j), 1e-2);
        }
    }
}
