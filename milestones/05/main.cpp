#include "verlet.h"                  // Includes the Verlet integration methods
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
    // Define a list of atom counts to test with the simulation
    std::vector<int> atom_counts = {50, 100, 200, 400, 800}; 
    std::list<double> time_store; // List to store simulation times for each atom count

    // Loop over different numbers of atoms
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

        // Set target temperature for the Berendsen thermostat
        double target_temperature = 100.0;
        double relaxation_time = 0.01;  // Relaxation time for temperature control

        // Open output files for trajectory data and energy values
        std::ofstream traj("../milestones/05/traj.xyz");             // File for trajectory
        std::ofstream energy_output("../milestones/05/energy.txt");  // File for energy values

        // Start timer to measure simulation time for this atom count
        auto start_time = std::chrono::high_resolution_clock::now();

        // Main simulation loop over time steps
        for (double step = 0; step < total_time; step += dt) {
            
            // Step 1 of Verlet integration: update positions and half-step velocities
            verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt);

            // Calculate potential energy using Lennard-Jones direct summation
            double potential_energy = lj_direct_summation(atoms);
            
            // Calculate kinetic energy: KE = 0.5 * sum(v^2)
            double kinetic_energy = 0.5 * atoms.velocities.square().sum();
            
            // Calculate total energy: E_total = KE + PE
            double total_energy = potential_energy + kinetic_energy;
            
            // Step 2 of Verlet integration: update velocities for the second half-step
            verlet_step2(atoms.velocities, atoms.forces, dt);

            // Apply the Berendsen thermostat to regulate temperature
            berendsen_thermostat(atoms, target_temperature, dt, relaxation_time);

            // Write the current atomic positions to the trajectory file
            write_xyz(traj, atoms);

            // Output energy information every 500 steps
            if (static_cast<int>(step) % 500 == 0) {
                std::cout << "Potential Energy: " << potential_energy << std::endl;
                std::cout << "Kinetic Energy: " << kinetic_energy << std::endl;
                std::cout << "Step: " << step << ", Total Energy: " << total_energy << std::endl;
            }

            // Write the current step and total energy to the energy output file
            energy_output << step * dt << " " << total_energy << std::endl;
        }

        // Stop the timer and calculate simulation time for this atom count
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> simulation_time = end_time - start_time;
        
        // Output the simulation time for the current atom count
        std::cout << "Simulation time for " << num_atoms << " atoms: " << simulation_time.count() << " seconds." << std::endl;

        // Store the simulation time
        time_store.push_back(simulation_time.count());

        // Close the trajectory and energy output files
        traj.close();
        energy_output.close();
    }

    // Output all the recorded simulation times
    std::cout << "Time store: ";
    for (const auto &time : time_store) {
        std::cout << time << " ";
    }
    std::cout << std::endl;

    return 0; // Program completed successfully
}
