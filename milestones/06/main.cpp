#include "verlet.h"
#include "lj.h"
#include "lj_direct_summation.h"
#include "xyz.h"
#include "berendsen_thermostat.h"
#include <list>
#include <fstream>
#include <iostream>
#include <cmath>
#include <chrono>  // Include chrono library for measuring time
#include <vector>  // Include vector for storing atom counts

// Function to create a simple cubic lattice of atoms
Positions_t create_cubic_lattice(int num_atoms, double lattice_constant) {
    int num_per_side = std::ceil(std::cbrt(num_atoms));
    Positions_t positions(3, num_atoms);
    int index = 0;
    for (int x = 0; x < num_per_side; ++x) {
        for (int y = 0; y < num_per_side; ++y) {
            for (int z = 0; z < num_per_side; ++z) {
                if (index < num_atoms) {
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
    std::cout<<"Cutoff: 5 used"<<std::endl;
    std::vector<int> atom_counts = {100,200,500,1000,1500,2000};
    std::list<double> time_store; // Store the simulation times for each atom count

    for (int num_atoms : atom_counts) {
        double lattice_constant = 1.2;
        Positions_t positions = create_cubic_lattice(num_atoms, lattice_constant);
        Velocities_t velocities(3, num_atoms);
        velocities.setZero();  // Initialize with zero velocities
        Atoms atoms(positions, velocities);

        double dt = 0.001 * sqrt(1 * pow(1, 2) / 1); // 0.001 * sqrt(m * sigma^2 / epsilon)
        double total_time = 100 * sqrt(1 * pow(1, 2) / 1); // 100 * sqrt(m * sigma^2 / epsilon)  // Adjust as necessary
        double target_temperature = 100.0;
        double relaxation_time =0.01;
        double cutoff = 5;


        auto start_time = std::chrono::high_resolution_clock::now();
        //for (double step = 0; step < total_time; step+=dt) {
        for (double step = 0; step < 100; step++) {


            verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt);
            double potential_energy = lj_cutoff_direct_summation(atoms, cutoff);
            verlet_step2(atoms.velocities, atoms.forces, dt);

            // Apply the Berendsen thermostat
            berendsen_thermostat(atoms, target_temperature, dt, relaxation_time);

            //write_xyz(traj, atoms);

            if (static_cast<int>(step) % 500 == 0) { // Output energy every 500 steps
               std::cout << "Potential Energy: " << potential_energy << std::endl;
           //     std::cout << "Kinetic Energy: " << kinetic_energy << std::endl;
           //     std::cout << "Step: " << step << ", Total Energy: " << total_energy << std::endl;
            }

        }

        auto end_time = std::chrono::high_resolution_clock::now();  // Stop measuring time
        std::chrono::duration<double> simulation_time = end_time - start_time;
        std::cout << "Simulation time for " << num_atoms << " atoms: " << simulation_time.count() << " seconds." << std::endl;
        time_store.push_back(simulation_time.count());

    }
    std::cout << "Time store: "; // Output the simulation times for each atom count
    for (const auto &time : time_store) {
        std::cout << time << " ";
    }
    std::cout << std::endl;

    std::cout<<"no cutoff used"<<std::endl;// Store the simulation times for each atom count
    time_store = {}; // Clear the time store
    for (int num_atoms : atom_counts) {
        double lattice_constant = 1.2;
        Positions_t positions = create_cubic_lattice(num_atoms, lattice_constant);
        Velocities_t velocities(3, num_atoms);
        velocities.setZero();  // Initialize with zero velocities
        Atoms atoms(positions, velocities);

        double dt = 0.001 * sqrt(1 * pow(1, 2) / 1); // 0.001 * sqrt(m * sigma^2 / epsilon)
        double total_time = 100 * sqrt(1 * pow(1, 2) / 1); // 100 * sqrt(m * sigma^2 / epsilon)  // Adjust as necessary
        double target_temperature = 100.0;
        double relaxation_time =0.01;


        auto start_time = std::chrono::high_resolution_clock::now();
        //for (double step = 0; step < total_time; step+=dt) { // Loop over time steps but as it was very time consuming had to decrease the timestep to 5000
        for (double step = 0; step < 5000; step++) {


            verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt);
            double potential_energy = lj_direct_summation(atoms);
            //std::cout << "Potential Energy: " << potential_energy << std::endl;
            double kinetic_energy = 0.5 * atoms.velocities.square().sum();
            //std::cout << "Kinetic Energy: " << kinetic_energy << std::endl;
            double total_energy = potential_energy + kinetic_energy;
            verlet_step2(atoms.velocities, atoms.forces, dt);

            // Apply the Berendsen thermostat
            berendsen_thermostat(atoms, target_temperature, dt, relaxation_time);



            if (static_cast<int>(step) % 500 == 0) { // Output energy every 500 steps
                std::cout << "Potential Energy: " << potential_energy << std::endl;
                std::cout << "Kinetic Energy: " << kinetic_energy << std::endl;
               std::cout << "Step: " << step << ", Total Energy: " << total_energy << std::endl;
            }

        }

        auto end_time = std::chrono::high_resolution_clock::now();  // Stop measuring time
        std::chrono::duration<double> simulation_time = end_time - start_time;
        std::cout << "Simulation time for " << num_atoms << " atoms: " << simulation_time.count() << " seconds." << std::endl;
        time_store.push_back(simulation_time.count());
    }
    std::cout << "Time store: "; // Output the simulation times for each atom count
    for (const auto &time : time_store) {
        std::cout << time << " ";
    }
    std::cout << std::endl;
    return 0;
}
