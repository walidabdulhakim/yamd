#include "verlet.h"
#include "lj_direct_summation.h"
#include "xyz.h"
#include <fstream>
#include <iostream>

int main() {
    auto [names, positions, velocities] = read_xyz_with_velocities("../milestones/04/lj54.xyz");
    Atoms atoms(positions, velocities);

    double dt = 0.001 * sqrt(1 * pow(1, 2) / 1); // 0.001 * sqrt(m * sigma^2 / epsilon)
    double total_time = 100 * sqrt(1 * pow(1, 2) / 1); // 100 * sqrt(m * sigma^2 / epsilon)  // Adjust as necessary

    std::ofstream traj("../milestones/04/traj.xyz");
    std::ofstream energy_output("../milestones/04/energy04.txt");
    double tolerance = 0.0001; // Adjust as necessary

    for (double step = 0; step < total_time; step +=dt ) {


        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt);
        double potential_energy = lj_direct_summation(atoms);
        double kinetic_energy = 0.5 * atoms.velocities.square().sum();
        double total_energy = potential_energy + kinetic_energy;
        verlet_step2(atoms.velocities, atoms.forces, dt);

        write_xyz(traj, atoms);

        if (std::abs(step - std::round(step)) < tolerance) {
            std::cout << "Step: " << step << ", Total Energy: " << total_energy << std::endl;
        }
        energy_output << step * dt << " " << total_energy << std::endl;
    }

    traj.close();
    energy_output.close();
    return 0;
}
