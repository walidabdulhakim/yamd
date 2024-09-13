#include "verlet.h"                  // Includes the Verlet integration methods
#include "lj.h"                      // Includes Lennard-Jones potential functions (possibly with a cutoff)
#include "lj_direct_summation.h"      // Includes Lennard-Jones potential calculation using direct summation
#include "xyz.h"                      // Includes methods for reading/writing .xyz files (for molecular data)
#include "berendsen_thermostat.h"     // Includes the Berendsen thermostat method
#include <list>                       // For storing simulation times
#include <fstream>                    // For file input/output operations
#include <iostream>                   // For standard I/O operations
#include <cmath>                      // For mathematical operations (e.g., sqrt, pow)
#include <chrono>                     // For timing the simulation
#include <vector>                     // For storing the list of atom counts

// Function to create a simple cubic lattice of atoms
// Creates atom positions arranged in a cubic grid with the given lattice constant
Positions_t create_cubic_lattice(int num_atoms, double lattice_constant) {
    int num_per_side = std::ceil(std::cbrt(num_atoms)); // Calculate the number of atoms per side of the cube
    Positions_t positions(3, num_atoms);                // Initialize position matrix (3D positions for each atom)
    int index = 0;
    // Loop to position atoms in a cubic grid
    for (int x = 0; x < num_per_side; ++x) {
        for (int y = 0; y < num_per_side; ++y) {
            for (int z = 0; z < num_per_side; ++z) {
                if (index < num_atoms) {
                    // Set atom positions on the cubic lattice, scaled by lattice_constant
                    positions(0, index) = x * lattice_constant;
                    positions(1, index) = y * lattice_constant;
                    positions(2, index) = z * lattice_constant;
                    ++index;
                }
            }
        }
    }
    return positions;
}

int main() {
    std::cout << "Cutoff: 5 used" << std::endl;  // Notify that a cutoff of 5 is being used
    std::vector<int> atom_counts = {100, 200, 500, 1000, 1500, 2000}; // List of atom counts to simulate
    std::list<double> time_store;  // Store the simulation times for each atom count

    // First simulation loop with Lennard-Jones potential using a cutoff
    for (int num_atoms : atom_counts) {
        double lattice_constant = 1.2; // Lattice constant for spacing atoms in the cubic grid
        Positions_t positions = create_cubic_lattice(num_atoms, lattice_constant); // Create cubic lattice positions
        Velocities_t velocities(3, num_atoms); // Initialize velocities (3D) for each atom
        velocities.setZero();  // Set initial velocities to zero (starting from rest)
        Atoms atoms(positions, velocities); // Create an Atoms object

        // Define time step (dt) using the formula: dt = 0.001 * sqrt(m * sigma^2 / epsilon)
        double dt = 0.001 * sqrt(1 * pow(1, 2) / 1); 

        // Define total simulation time: total_time = 100 * sqrt(m * sigma^2 / epsilon)
        double total_time = 100 * sqrt(1 * pow(1, 2) / 1); 

        // Parameters for Berendsen thermostat
        double target_temperature = 100.0;
        double relaxation_time = 0.01;  // Relaxation time for temperature control
        double cutoff = 5.0;  // Lennard-Jones potential cutoff distance

        // Start timer to measure simulation time for this atom count
        auto start_time = std::chrono::high_resolution_clock::now();

        // Main simulation loop with Lennard-Jones cutoff potential
        for (double step = 0; step < 100; step++) { // Loop for a fixed number of steps (100 steps)

            // Step 1 of Verlet integration: update positions and half-step velocities
            verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt);

            // Calculate potential energy using Lennard-Jones direct summation with a cutoff
            double potential_energy = lj_cutoff_direct_summation(atoms, cutoff);

            // Step 2 of Verlet integration: update velocities for the second half-step
            verlet_step2(atoms.velocities, atoms.forces, dt);

            // Apply the Berendsen thermostat to regulate temperature
            berendsen_thermostat(atoms, target_temperature, dt, relaxation_time);

            // Output energy information every 500 steps (if the simulation were larger)
            if (static_cast<int>(step) % 500 == 0) {
                std::cout << "Potential Energy: " << potential_energy << std::endl;
            }
        }

        // Stop the timer and calculate simulation time for this atom count
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> simulation_time = end_time - start_time;
        
        // Output the simulation time for the current atom count
        std::cout << "Simulation time for " << num_atoms << " atoms: " << simulation_time.count() << " seconds." << std::endl;

        // Store the simulation time
        time_store.push_back(simulation_time.count());
    }

    // Output all the recorded simulation times for the cutoff run
    std::cout << "Time store: ";
    for (const auto &time : time_store) {
        std::cout << time << " ";
    }
    std::cout << std::endl;

    // Second simulation without Lennard-Jones cutoff
    std::cout << "no cutoff used" << std::endl;  // Notify that no cutoff is being used
    time_store = {};  // Clear the time store to reuse for the second set of simulations

    // Repeat the same simulation but without a cutoff for the Lennard-Jones potential
    for (int num_atoms : atom_counts) {
        double lattice_constant = 1.2;
        Positions_t positions = create_cubic_lattice(num_atoms, lattice_constant);
        Velocities_t velocities(3, num_atoms);
        velocities.setZero();
        Atoms atoms(positions, velocities);

        double dt = 0.001 * sqrt(1 * pow(1, 2) / 1);
        double total_time = 100 * sqrt(1 * pow(1, 2) / 1);
        double target_temperature = 100.0;
        double relaxation_time = 0.01;

        auto start_time = std::chrono::high_resolution_clock::now();

        // Main simulation loop without cutoff
        for (double step = 0; step < 5000; step++) { // Increase the number of steps to 5000

            // Step 1 of Verlet integration
            verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt);

            // Calculate potential energy using Lennard-Jones direct summation without cutoff
            double potential_energy = lj_direct_summation(atoms);

            // Calculate kinetic energy
            double kinetic_energy = 0.5 * atoms.velocities.square().sum();

            // Calculate total energy (potential + kinetic)
            double total_energy = potential_energy + kinetic_energy;

            // Step 2 of Verlet integration
            verlet_step2(atoms.velocities, atoms.forces, dt);

            // Apply the Berendsen thermostat
            berendsen_thermostat(atoms, target_temperature, dt, relaxation_time);

            // Output energy information every 500 steps
            if (static_cast<int>(step) % 500 == 0) {
                std::cout << "Potential Energy: " << potential_energy << std::endl;
                std::cout << "Kinetic Energy: " << kinetic_energy << std::endl;
                std::cout << "Step: " << step << ", Total Energy: " << total_energy << std::endl;
            }
        }

        // Stop the timer and calculate simulation time
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> simulation_time = end_time - start_time;
        std::cout << "Simulation time for " << num_atoms << " atoms: " << simulation_time.count() << " seconds." << std::endl;

        // Store the simulation time
        time_store.push_back(simulation_time.count());
    }

    // Output all the recorded simulation times for the no-cutoff run
    std::cout << "Time store: ";
    for (const auto &time : time_store) {
        std::cout << time << " ";
    }
    std::cout << std::endl;

    return 0;  // Program completed successfully
}
