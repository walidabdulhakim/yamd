#include "verlet.h"                 // Includes the Verlet integration methods
#include "lj_direct_summation.h"     // Includes the Lennard-Jones direct summation potential calculation methods
#include "xyz.h"                     // Includes methods for reading/writing .xyz files (for molecular data)
#include <fstream>                   // For file input/output operations
#include <iostream>                  // For standard I/O operations

int main() {
    // Read atom names, positions, and velocities from the xyz file
    auto [names, positions, velocities] = read_xyz_with_velocities("../milestones/04/lj54.xyz");
    
    // Create an Atoms object to store atomic positions and velocities
    Atoms atoms(positions, velocities);

    // Define time step (dt) using the formula: dt = 0.001 * sqrt(m * sigma^2 / epsilon)
    // For this simulation, m = sigma = epsilon = 1, so dt simplifies to 0.001
    double dt = 0.001 * sqrt(1 * pow(1, 2) / 1); 

    // Define total simulation time using the formula: total_time = 100 * sqrt(m * sigma^2 / epsilon)
    double total_time = 100 * sqrt(1 * pow(1, 2) / 1); 

    // Open output files for trajectory data (in .xyz format) and energy values
    std::ofstream traj("../milestones/04/traj.xyz");             // File for trajectory
    std::ofstream energy_output("../milestones/04/energy04.txt"); // File for energy values

    // Define tolerance for printing step information (used for checking when to output step info)
    double tolerance = 0.0001; // Can be adjusted as necessary

    // Main loop for time integration using Verlet algorithm
    for (double step = 0; step < total_time; step += dt) {
        
        // Step 1 of Verlet integration: update positions and half-step velocities
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt);

        // Calculate potential energy using Lennard-Jones direct summation method
        double potential_energy = lj_direct_summation(atoms);

        // Calculate kinetic energy: KE = 0.5 * sum(v^2)
        double kinetic_energy = 0.5 * atoms.velocities.square().sum();

        // Calculate total energy: E_total = KE + PE
        double total_energy = potential_energy + kinetic_energy;

        // Step 2 of Verlet integration: update velocities for the second half-step
        verlet_step2(atoms.velocities, atoms.forces, dt);

        // Write the current atomic positions to the trajectory file
        write_xyz(traj, atoms);

        // Output step information to the console if the step is close to an integer value
        if (std::abs(step - std::round(step)) < tolerance) {
            std::cout << "Step: " << step << ", Total Energy: " << total_energy << std::endl;
        }

        // Write the current step time and total energy to the energy output file
        energy_output << step * dt << " " << total_energy << std::endl;
    }

    // Close the trajectory and energy output files
    traj.close();
    energy_output.close();

    return 0; // Program completed successfully
}
