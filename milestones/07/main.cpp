#include "atoms.h"
#include "verlet.h"
#include "xyz.h"
#include <ducastelle.h>
#include <iomanip>
#include <iostream>
#include <neighbors.h>
#include <berendsen_thermostat.h>

using namespace std;
// cutoff distance of the ducastelle potential
const double CUTOFF_DISTANCE{10.};

const vector<string> CL_files{

    "../mackay_gold/cluster_6.xyz",
    "../mackay_gold/cluster_7.xyz",
    "../mackay_gold/cluster_8.xyz",
    "../mackay_gold/cluster_9.xyz",
    "../mackay_gold/cluster_10.xyz",

};

/*A helper function to subtract the mean velocity for each velocity component.
 * This reduces code duplication and enhances readability.*/
void zero_net_velocity(Atoms &atoms) {
    for (int i = 0; i < atoms.velocities.rows(); ++i) {
        atoms.velocities.row(i) -= atoms.velocities.row(i).mean();
    }
}

/*heats up the atom to desired target temp*/
void heat_up(Atoms &atoms, double dt, double target_temp) {
    // Initialize atoms with a small random initial velocity
    atoms.velocities.setRandom();
    atoms.velocities *= 1e-3;

    // Create neighbor list and initialize forces
    NeighborList neighbor_list;
    neighbor_list.update(atoms, CUTOFF_DISTANCE);
    if (!ducastelle(atoms, neighbor_list)) {
        cerr << "Error in ducastelle force calculation" << endl;
        return;
    }

    double time_elapsed = 0;
    const double relaxation_time = 1000 * dt;
    const double heating_duration = 20 * relaxation_time;
    const double rest_duration = 10 * relaxation_time;

    while (time_elapsed < heating_duration + rest_duration) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt, atoms.masses);
        neighbor_list.update(atoms, CUTOFF_DISTANCE);
        if (!ducastelle(atoms, neighbor_list)) {
            cerr << "Error in ducastelle force calculation" << endl;
            return;
        }
        verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);

        if (time_elapsed < heating_duration) {
            berendsen_thermostat(atoms, target_temp, dt, relaxation_time );
            zero_net_velocity(atoms);
        }

        cout << "T=" << temperature_cur(atoms) << "K\r" << flush;

        time_elapsed += dt;
    }

    cout << "T=" << temperature_cur(atoms) << "K" << endl;
}


double get_potenergy(Atoms &atoms, double dt, NeighborList &neighbour_list) {
    verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt,
                 atoms.masses);
    neighbour_list.update(atoms, CUTOFF_DISTANCE);
    double potential = ducastelle(atoms, neighbour_list);
    //cout << "E_pot=" << potential <<  endl;
    verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);
    return potential;
};



int main(int argc, char *argv[]) {

    if (argc > 1 && string(argv[1]) == "H") {
        auto [names, positions,
              velocities]= read_xyz_with_velocities("../milestones/07/cluster_923.xyz");
        Atoms atoms= Atoms(names, positions);
        // heat up the system to 500K and save the result to a file that's
        heat_up(atoms, 0.1, 500.);
        ofstream traj("../milestones/07/923_500K.xyz");
        write_xyz(traj, atoms);
        return 0;
    }

    // find the optimal timestep size (.1 is a very stable one for this system
    else if (argc > 1 && string(argv[1]) == "T") {
        // total time to simulate the system for
        // double dt = 0.001 * sqrt(1 * pow(1, 2) / 1); // 0.001 * sqrt(m * sigma^2 / epsilon)
        double total_time = 100 * sqrt(1 * pow(1, 2) / 1);
        // open a csv file for the results, write the header
        ofstream dt_eam_file("../milestones/07/dt_eam.csv");
        dt_eam_file << "time step (dt), time (t) ,Total Energy, potential energy(E_pot) , kinetic energy (E_kin),"
                    << endl;

        for (double dt : {10., 1., 0.5, 0.1}) {
            const string filename = "../milestones/07/923_500K.xyz";
            auto [names, positions, velocities] =
                read_xyz_with_velocities(filename);
            Atoms atoms = Atoms(names, positions);

            // create neighbour list and initialize forces
            NeighborList neighbor_list;
            neighbor_list.update(atoms, CUTOFF_DISTANCE);
            ducastelle(atoms, neighbor_list);
            double t = 0;
            while (t < total_time) {
                double potential = get_potenergy(atoms, dt, neighbor_list);
                double kinetic = atoms.kinetic_energy();
                cout << "ke=" << kinetic << endl;
                dt_eam_file << dt << "," << t << "," << kinetic + potential
                            << "," << potential << "," << kinetic << endl;
                t += dt;
            }

        }
    }


        // EAM ISOCAHEDRON Energy temprature readings

    double dt{.5};// stable timestep
    size_t tau_relax{1000}; // in units of dt
    double from_temp{500};// initial temperature
    double until_temp{1000};
    double increase_temp_by{10.};


        // initialize csv file for the results
    ofstream gold_r("../milestones/07/Readings.csv");
    gold_r << "Cluster size,total energy,temperature,potential energy,kinetic energy"
           << endl;

        // loop over all clusters, heat them up and record the results
    for (string filename : CL_files) {
        // initialize gold cluster from file
         auto [names, positions,velocities]= read_xyz_with_velocities(filename);
         Atoms atoms= Atoms(names, positions);

         //ofstream traj("/mnt/c/Users/walid/OneDrive/Desktop/yamd/milestones/07/923_" +
                       //to_string(atoms.nb_atoms()) + ".xyz");

            //  pre-het the gold to some initial temperature
         cout << "Heating up cluster of size " << atoms.nb_atoms()
              << " to " << from_temp << "K" << endl;
         heat_up(atoms, dt, from_temp);
         cout << "Start of Reading the results" << endl;

         // create neighbour list and initialize forces
         NeighborList neighbour_list;
         neighbour_list.update(atoms, CUTOFF_DISTANCE);
         ducastelle(atoms, neighbour_list);
         double current_temp= 0.;
         double e_pot= 0.0;

         // // might as well record the trajectory to check plausibility of
         // results std::ofstream traj("gold_t_e_traj.xyz");

            // then, repeatedly increase temperature, wait for tau_relax, measure energy and temperature for tau_relax
         while (current_temp < until_temp) {
                // bump up energy
                berendsen_thermostat(atoms, temperature_cur(atoms) + increase_temp_by,1.0,1.0); // increase temperature by 10K
                cout << "temperature" << temperature_cur(atoms) << "K" << endl;
                cout<<"kinetic energy"<<atoms.kinetic_energy()<<endl;
                // wait for system to relax
                size_t counter = tau_relax;
                for (size_t i{0}; i <= tau_relax; i++)  {
                    verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt,
                                 atoms.masses);
                    neighbour_list.update(atoms, CUTOFF_DISTANCE);
                    double pot = ducastelle(atoms, neighbour_list);
                    //cout << "E_pot=" << pot << endl;
                    verlet_step2(atoms.velocities, atoms.forces, dt, atoms.masses);
                }
                //cout << "temperature" << temperature_cur(atoms) << "K" << endl;
                //cout<<"kinetic energy"<<atoms.kinetic_energy()<<endl;
                // measure T and E
                double avg_e_tot{0.0};
                double avg_e_pot{0.0};
                double avg_e_kin{0.0};
                double avg_t{0.0};
                for (size_t i{0}; i <= tau_relax; i++) {
                    e_pot = get_potenergy(atoms, dt, neighbour_list);
                    double e_kin = atoms.kinetic_energy();
                    avg_e_tot += e_pot + e_kin;
                    avg_e_pot += e_pot;
                    avg_e_kin += e_kin;
                    avg_t += temperature_cur(atoms);
                }
                // write data to csv
                double normalize = 1. / ((double)tau_relax);
                gold_r << atoms.nb_atoms() << "," << fixed
                         << setprecision(15) << avg_e_tot * normalize
                         << "," << avg_t * normalize << ","
                         << avg_e_pot * normalize << ","
                         << avg_e_kin * normalize << endl;

                //write_xyz(traj,atoms);// write trajectory to file for later inspection

                // update current temperature
                current_temp = temperature_cur(atoms);


                cout << atoms.nb_atoms()
                          << " atoms heated up to: " << avg_t * normalize
                          << "K | E_total: " << avg_e_tot * normalize
                          << "eV | E_pot: " << avg_e_pot * normalize
                          << "eV | E_kin: " << avg_e_kin * normalize << "eV"
                          << endl;
            }
        }

    return 0;
}
