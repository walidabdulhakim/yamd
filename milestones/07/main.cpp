#include "atoms.h"                    // Defines the Atoms class and related functionalities
#include "verlet.h"                   // Includes the Verlet integration methods
#include "xyz.h"                      // Includes methods for reading/writing .xyz files (for molecular data)
#include <ducastelle.h>               // Includes the Ducastelle potential for energy and force calculations
#include <iomanip>                    // For controlling output formatting
#include <iostream>                   // For standard I/O operations
#include <neighbors.h>                // For neighbor list generation (used in potentials)
#include <berendsen_thermostat.h>     // Includes the Berendsen thermostat method

using namespace std;

// Define the cutoff distance for the Ducastelle potential
const double CUTOFF_DISTANCE{10.};

// List of .xyz files for different cluster sizes
const vector<string> CL_files{
    "../mackay_gold/cluster_6.xyz",
    "../mackay_gold/cluster_7.xyz",
    "../mackay_gold/cluster_8.xyz",
    "../mackay_gold/cluster_9.xyz",
    "../mackay_gold/cluster_10.xyz",
};

/* Function to zero the net velocity of the system by subtracting the mean velocity.
 * This ensures that the system's center of mass is not moving. */
void zero_net_velocity(Atoms &atoms) {
    for (int i = 0; i < atoms.velocities.rows(); ++i) {
        atoms.velocities.row(i) -= atoms.velocities.row(i).mean();
    }
}

/* Function to heat the atoms to a target temperature using the Berendsen thermostat. */
void heat_up(Atoms &atoms, double dt, double target_temp) {
    // Initialize atom velocities with small random values
    atoms.velocities.setRandom();
    atoms.velocities *= 1e-3;

    // Create a neighbor list for efficient force calculation
    NeighborList neighbor_list;
    neighbor_list.update(atoms, CUTOFF_DISTANCE);

    // Calculate initial forces using the Ducastelle potential
    if (!ducastelle(atoms, neighbor_list)) {
        cerr << "Error in ducastelle force calculation" << endl;
        return;
    }

    // Time-related parameters for the heating process
    double time_elapsed = 0;
    const double relaxation_time = 1000 * dt;
    const double heating_duration = 20 * relaxation_time;
    const double rest_duration = 10 * relaxation_time;

    // Loop over time to simulate the heating and relaxation
    while (time_elapsed < heating_duration + rest_duration) {
        // Perform one Verlet integration step
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt, atoms.masses);
        neighbor_list.update(atoms, CUTOFF_DISTANCE);
        if (!ducastelle(atoms, neighbor_list)) {
            cerr << "Error in ducastelle force calculation" << endl;
            return;
        }
        verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);

        // Apply the Berendsen thermostat during the heating phase
        if (time_elapsed < heating_duration) {
            berendsen_thermostat(atoms, target_temp, dt, relaxation_time);
            zero_net_velocity(atoms);  // Ensure the system's center of mass isn't moving
        }

        // Output the current temperature
        cout << "T=" << temperature_cur(atoms) << "K\r" << flush;
        time_elapsed += dt;
    }

    // Final temperature after heating
    cout << "T=" << temperature_cur(atoms) << "K" << endl;
}

/* Function to perform a single step of Verlet integration and return the potential energy */
double get_potenergy(Atoms &atoms, double dt, NeighborList &neighbor_list) {
    verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt, atoms.masses);
    neighbor_list.update(atoms, CUTOFF_DISTANCE);
    double potential = ducastelle(atoms, neighbor_list);  // Calculate potential energy using Ducastelle potential
    verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);
    return potential;
}

int main(int argc, char *argv[]) {
    // Case for heating the system to 500K and saving the result to a file
    if (argc > 1 && string(argv[1]) == "H") {
        auto [names, positions, velocities] = read_xyz_with_velocities("../milestones/07/cluster_923.xyz");
        Atoms atoms = Atoms(names, positions);

        // Heat the system to 500K
        heat_up(atoms, 0.1, 500.);

        // Save the heated configuration to an XYZ file
        ofstream traj("../milestones/07/923_500K.xyz");
        write_xyz(traj, atoms);
        return 0;
    }

    // Case for testing different timestep sizes and recording energy values
    else if (argc > 1 && string(argv[1]) == "T") {
        double total_time = 100 * sqrt(1 * pow(1, 2) / 1);  // Define total simulation time

        // Open a CSV file to record the results
        ofstream dt_eam_file("../milestones/07/dt_eam.csv");
        dt_eam_file << "time step (dt), time (t), Total Energy, potential energy (E_pot), kinetic energy (E_kin)" << endl;

        // Loop over different time steps
        for (double dt : {10., 1., 0.5, 0.1}) {
            const string filename = "../milestones/07/923_500K.xyz";
            auto [names, positions, velocities] = read_xyz_with_velocities(filename);
            Atoms atoms = Atoms(names, positions);

            // Create neighbor list and initialize forces
            NeighborList neighbor_list;
            neighbor_list.update(atoms, CUTOFF_DISTANCE);
            ducastelle(atoms, neighbor_list);
            double t = 0;

            // Simulate the system for the specified time
            while (t < total_time) {
                double potential = get_potenergy(atoms, dt, neighbor_list);
                double kinetic = atoms.kinetic_energy();

                // Record energy values in the CSV file
                dt_eam_file << dt << "," << t << "," << kinetic + potential << "," << potential << "," << kinetic << endl;
                t += dt;
            }
        }
    }

    // Parameters for heating the gold clusters
    double dt{0.5};            // Stable timestep
    size_t tau_relax{1000};     // Relaxation time in units of dt
    double from_temp{500};      // Initial temperature (500K)
    double until_temp{1000};    // Target temperature (1000K)
    double increase_temp_by{10.}; // Temperature increment (10K)

    // Initialize CSV file to store the results of cluster heating
    ofstream gold_r("../milestones/07/Readings.csv");
    gold_r << "Cluster size, total energy, temperature, potential energy, kinetic energy" << endl;

    // Loop over all cluster files to heat them up and record results
    for (string filename : CL_files) {
        // Load the atom data from the file
        auto [names, positions, velocities] = read_xyz_with_velocities(filename);
        Atoms atoms = Atoms(names, positions);

        // Heat the cluster to the initial temperature (500K)
        cout << "Heating up cluster of size " << atoms.nb_atoms() << " to " << from_temp << "K" << endl;
        heat_up(atoms, dt, from_temp);
        cout << "Start of Reading the results" << endl;

        // Create neighbor list and initialize forces
        NeighborList neighbor_list;
        neighbor_list.update(atoms, CUTOFF_DISTANCE);
        ducastelle(atoms, neighbor_list);
        double current_temp = 0.0;
        double e_pot = 0.0;

        // Increment the temperature and relax the system at each step
        while (current_temp < until_temp) {
            // Increase temperature by 10K using the Berendsen thermostat
            berendsen_thermostat(atoms, temperature_cur(atoms) + increase_temp_by, 1.0, 1.0);
            cout << "Temperature: " << temperature_cur(atoms) << "K" << endl;
            cout << "Kinetic energy: " << atoms.kinetic_energy() << endl;

            // Wait for the system to relax
            for (size_t i = 0; i <= tau_relax; i++) {
                verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt, atoms.masses);
                neighbor_list.update(atoms, CUTOFF_DISTANCE);
                double pot = ducastelle(atoms, neighbor_list);
                verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);
            }

            // Calculate and record the average temperature and energy
            double avg_e_tot = 0.0, avg_e_pot = 0.0, avg_e_kin = 0.0, avg_t = 0.0;
            for (size_t i = 0; i <= tau_relax; i++) {
                e_pot = get_potenergy(atoms, dt, neighbor_list);
                double e_kin = atoms.kinetic_energy();
                avg_e_tot += e_pot + e_kin;
                avg_e_pot += e_pot;
                avg_e_kin += e_kin;
                avg_t += temperature_cur(atoms);
            }

            // Normalize and write the results to the CSV file
            double normalize = 1.0 / ((double)tau_relax);
            gold_r << atoms.nb_atoms() << "," << fixed << setprecision(15) << avg_e_tot * normalize << ","
                   << avg_t * normalize << "," << avg_e_pot * normalize << "," << avg_e_kin * normalize << endl;

            // Update the current temperature for the next iteration
            current_temp = temperature_cur(atoms);

            // Output the results to the console
            cout << atoms.nb_atoms() << " atoms heated up to: " << avg_t * normalize
                 << "K | E_total: " << avg_e_tot * normalize << "eV | E_pot: " << avg_e_pot * normalize
                 << "eV | E_kin: " << avg_e_kin * normalize << "eV" << endl;
        }
    }

    return 0;  // Program completed successfully
}
