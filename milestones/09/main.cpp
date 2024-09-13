#include "atoms.h"                    // Defines the Atoms class and related functionalities
#include "mpi.h"                      // MPI (Message Passing Interface) for parallel computing
#include "mpi_support.h"              // Includes MPI support functions
#include "verlet.h"                   // Includes the Verlet integration methods
#include "xyz.h"                      // Includes methods for reading/writing .xyz files
#include <domain.h>                   // Domain decomposition for parallel processing
#include <ducastelle.h>               // Ducastelle potential for energy and force calculations
#include <iomanip>                    // For controlling output formatting
#include <iostream>                   // For standard I/O operations
#include <neighbors.h>                // Neighbor list generation (used in potentials)
#include <berendsen_thermostat.h>     // Includes the Berendsen thermostat method

// Define constants for atomic spacing and border size
const double spacing = 4.079 / sqrt(2); // Spacing between atoms
const double border = 50 * spacing;     // Border size for volume calculations

/* Function to remove the mean velocity from the system by subtracting the average velocity.
 * Ensures the center of mass of the system is stationary. */
void remove_mean_velocity(Velocities_t &v) {
    Vec3_t mean_velocity = v.rowwise().mean();
    v.row(0) -= mean_velocity(0);
    v.row(1) -= mean_velocity(1);
    v.row(2) -= mean_velocity(2);
}

/* Function to calculate the volume of the system based on the atom positions.
 * It adjusts the positions and adds border spacing to get the dimensions. */
Vec3_t calculate_volume(const Atoms &atoms) {
    Vec3_t offset{atoms.positions.rowwise().minCoeff()}; // Find minimum position for offset

    // Adjust positions by subtracting the offset
    auto adjusted_positions = atoms.positions;
    adjusted_positions.row(0) -= offset(0);
    adjusted_positions.row(1) -= offset(1);
    adjusted_positions.row(2) -= offset(2);

    // Calculate the minimum and maximum positions after adjustment
    Vec3_t min_pos = adjusted_positions.rowwise().minCoeff();
    Vec3_t max_pos = adjusted_positions.rowwise().maxCoeff();

    // Compute the volume dimensions
    Vec3_t volume_dimensions = max_pos - min_pos;

    // Add border and spacing adjustments
    volume_dimensions(0) += 2.0 * border;
    volume_dimensions(1) += 2.0 * border;
    volume_dimensions(2) += spacing / 4.0;

    return volume_dimensions;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);  // Initialize MPI for parallel processing

    // Simulation parameters
    const double dt = 5;                          // Time step for simulation
    const double TAU = 1000.;                     // Relaxation time for the Berendsen thermostat
    const size_t Scale = (size_t)round(10 * TAU); // Scaling factor based on relaxation time
    const double temperature[] = {0., 600.};      // Temperature settings (e.g., 0K, 600K)
    const double strainRateHz[] = {1e8, 5e8};     // Strain rates in Hz
    const double stretchByPercent = 40;           // Percentage stretch applied to the system

    // Filename for input data
    const std::string FILENAME = "../milestones/09/whisker_small";

    // Loop over different temperatures and strain rates
    for (auto temp : temperature) {
        for (auto strainRate : strainRateHz) {
            std::cout << "temp: " << temp << " strainRate: " << strainRate << std::endl;

            // Simulation constants
            const double cutoff = 10.;  // Cutoff distance for the Ducastelle potential
            const size_t TIMESTEPS = (size_t)ceil((stretchByPercent / 100.) /
                                                   (strainRate) / (dt * 1e-15)); // Total simulation time steps
            const size_t plotEvery = TIMESTEPS / 1000;  // Plot every 1000 steps

            // Read atom names, positions, and velocities from XYZ file
            auto [names, positions, velocities] = read_xyz_with_velocities(FILENAME + ".xyz");
            Atoms atoms{Atoms(names, positions)};  // Initialize Atoms object
            atoms.velocities.setRandom();          // Set random initial velocities
            atoms.velocities *= 1e-3;              // Scale down the velocities
            remove_mean_velocity(atoms.velocities); // Ensure no net velocity

            // Mass of gold (Au) atoms
            double m = ELEM_NAME_TO_MASS.at("Au");

            // Open files for writing trajectory and results
            std::ofstream traj(FILENAME + "_" + std::to_string(temp) + "_" + std::to_string(strainRate) + "_traj.xyz");
            std::ofstream results(FILENAME + "_" + std::to_string(temp) + "_" + std::to_string(strainRate) + ".csv");

            // Initialize neighbor list and calculate initial forces
            NeighborList neighbour_list;
            neighbour_list.update(atoms, cutoff);
            ducastelle(atoms, neighbour_list);

            int nb_processes = 0;
            MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);  // Get the number of MPI processes

            // Calculate the volume of the system
            Vec3_t volume = calculate_volume(atoms);

            // Initialize domain decomposition for parallel computation
            Domain domain(MPI_COMM_WORLD, volume, {1, 1, nb_processes}, {0, 0, 1});

            // Write header to the results CSV file
            if (domain.rank() == 0) {
                results << "Temperature,force_zz,volume_z,Potential_energy,Time,v_z_0,strainRateHz,volume_x,volume_y" << std::endl;
            }

            // Enable domain decomposition for the atoms
            domain.enable(atoms);

            // Update ghost atoms for inter-process communication
            domain.update_ghosts(atoms, cutoff * 2.);
            neighbour_list.update(atoms, cutoff);
            ducastelle(atoms, neighbour_list);

            // Time-stepping loop
            double t = 0;                   // Current simulation time
            double force_local = 0.;        // Local force accumulation
            double v_z_0 = volume(2);       // Initial Z-dimension of the volume
            for (size_t i = 0; i < TIMESTEPS + Scale; i++) {
                if (i >= Scale) {
                    if (i == Scale) {
                        // Recalculate volume and enable domain scaling
                        domain.disable(atoms);
                        volume = calculate_volume(atoms);
                        v_z_0 = volume(2);
                        domain.enable(atoms);
                    }
                    // Apply strain to the Z-dimension of the volume
                    volume(2) = v_z_0 + v_z_0 * strainRate * t * (1e-15);
                    domain.scale(atoms, volume);  // Scale the domain according to the applied strain
                }

                // Perform Verlet integration steps
                verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt, m);
                domain.exchange_atoms(atoms);  // Exchange atom data between processes
                domain.update_ghosts(atoms, cutoff * 2.);
                neighbour_list.update(atoms, cutoff);

                // Calculate potential energy and stress
                auto [e_pot_l, stress_zz_l] = ducastelle_stress(atoms, neighbour_list, domain.nb_local(), cutoff);

                // Complete the second step of the Verlet integration
                verlet_step2(atoms.velocities, atoms.forces, dt, m);
                
                // Apply Berendsen thermostat to control temperature
                berendsen_thermostat_decomposed(atoms, temp, dt, TAU, domain.nb_local(), m);

                // Accumulate the local force in the Z direction
                if (i >= Scale) {
                    force_local += stress_zz_l(2, 2) / volume(2);
                }

                // Periodically output the results
                if (i >= Scale && (i - Scale) % plotEvery == 0) {
                    // Reduce and gather results from all processes
                    double force_zz_g = MPI::allreduce(force_local / plotEvery, MPI_SUM, domain.communicator());
                    force_local = 0.0;
                    double e_pot_global = MPI::allreduce(e_pot_l, MPI_SUM, domain.communicator());

                    domain.disable(atoms);   // Temporarily disable domain decomposition
                    remove_mean_velocity(atoms.velocities);  // Remove mean velocity

                    // Write results to the output files
                    if (domain.rank() == 0) {
                        std::cout << std::fixed << std::setprecision(15)
                                  << temperature_cur(atoms) << "," << force_zz_g << "," << volume(2)
                                  << "," << e_pot_global << "," << t << "," << v_z_0 << ","
                                  << strainRate << "," << volume(0) << "," << volume(1)
                                  << std::endl;

                        results << std::fixed << std::setprecision(15)
                                << temperature_cur(atoms) << "," << force_zz_g << "," << volume(2)
                                << "," << e_pot_global << "," << t << "," << v_z_0 << ","
                                << strainRate << "," << volume(0) << "," << volume(1)
                                << std::endl;

                        write_xyz(traj, atoms);  // Save trajectory data to the file
                    }

                    domain.enable(atoms);  // Re-enable domain decomposition
                    domain.update_ghosts(atoms, cutoff * 2.);
                    neighbour_list.update(atoms, cutoff);
                    ducastelle(atoms, neighbour_list);
                }

                if (i >= Scale) {
                    t += dt;  // Increment simulation time
                }
            }
        }
    }
    MPI_Finalize();  // Finalize MPI before program exit
    return 0;
}
