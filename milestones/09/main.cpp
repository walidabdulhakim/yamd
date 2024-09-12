#include "atoms.h"
#include "mpi.h"
#include "mpi_support.h"
#include "verlet.h"
#include "xyz.h"
#include <domain.h>
#include <ducastelle.h>
#include <iomanip>
#include <iostream>
#include <neighbors.h>
#include <berendsen_thermostat.h>



const double spacing =4.079 / sqrt(2);
const double border =50 * spacing;


void remove_mean_velocity(Velocities_t &v) {
    Vec3_t mean_velocity=v.rowwise().mean();
    v.row(0) -= mean_velocity(0);
    v.row(1) -= mean_velocity(1);
    v.row(2) -= mean_velocity(2);
}

Vec3_t calculate_volume(const Atoms &atoms) {
    // Calculate the minimum position
    Vec3_t offset{atoms.positions.rowwise().minCoeff()};

    // Adjust positions based on the minimum position
    auto adjusted_positions = atoms.positions;
    adjusted_positions.row(0) -= offset(0);
    adjusted_positions.row(1) -= offset(1);
    adjusted_positions.row(2) -= offset(2);

    // Recalculate the minimum position after adjustment
    Vec3_t min_pos=adjusted_positions.rowwise().minCoeff();
    // Calculate the maximum position
    Vec3_t max_pos=adjusted_positions.rowwise().maxCoeff();

    // Determine the volume vector
    Vec3_t volume_dimensions=max_pos - min_pos;

    // Apply the border and spacing adjustments
    volume_dimensions(0) += 2.0 * border;
    volume_dimensions(1) += 2.0 * border;
    volume_dimensions(2) += spacing / 4.0;

    return volume_dimensions;
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // simulation settings
    const double dt=5;
    const double TAU=1000.;
    const size_t Scale=(size_t)round(10 * TAU);
    const double temperature[]={0.,600.};
    const double strainRateHz[] = {1e8, 5e8};
    const double stretchByPercent=40;

    const std::string FILENAME="../milestones/09/whisker_small";

    for (auto temp : temperature) {
        for (auto strainRate : strainRateHz) {
        std::cout << "temp: " << temp << " strainRate: " << strainRate << std::endl;
        // other constants
        const double cutoff=10.;
        const size_t TIMESTEPS=(size_t)ceil((stretchByPercent / 100.) /
                                              (strainRate) / (dt * 1e-15));
        const size_t plotEvery =TIMESTEPS / 1000;


        auto [names, positions, velocities]{read_xyz_with_velocities(FILENAME + ".xyz")};
        Atoms atoms{Atoms(names, positions)};
        atoms.velocities.setRandom();
        atoms.velocities *= 1e-3;
        remove_mean_velocity(atoms.velocities);

        double m = ELEM_NAME_TO_MASS.at("Au");
        std::ofstream traj(FILENAME +"_"+ std::to_string(temp)+"_"+std::to_string(strainRate) +"_traj.xyz");
        std::ofstream results(FILENAME +"_"+ std::to_string(temp)+"_"+std::to_string(strainRate) + ".csv");

        NeighborList neighbour_list;
        neighbour_list.update(atoms, cutoff);
        ducastelle(atoms, neighbour_list);

        int nb_processes=0;
        MPI_Comm_size(MPI_COMM_WORLD, &nb_processes);

        Vec3_t volume=calculate_volume(atoms);

        Domain domain(MPI_COMM_WORLD, volume, {1, 1, nb_processes}, {0, 0, 1});

        if (domain.rank() == 0) {
            results << "Temperature,force_zz,volume_z,Potential_energy,Time,v_z_0,strainRateHz,volume_x,volume_y" << std::endl;
        }

        domain.enable(atoms);

        domain.update_ghosts(atoms, cutoff * 2.);
        neighbour_list.update(atoms, cutoff);
        ducastelle(atoms, neighbour_list);

        double t=0;
        double force_local=0.;
        double v_z_0 = volume(2);
        for (size_t i=0; i < TIMESTEPS + Scale; i++) {
            if (i >= Scale) {
                if (i == Scale) {
                    domain.disable(atoms);
                    volume = calculate_volume(atoms);
                    v_z_0 = volume(2);
                    domain.enable(atoms);
                }
                volume(2) = v_z_0 + v_z_0 * strainRate * t * (1e-15);
                domain.scale(atoms, volume);
            }
            verlet_step1(atoms.positions, atoms.velocities, atoms.forces, dt, m);
            domain.exchange_atoms(atoms);
            domain.update_ghosts(atoms, cutoff * 2.);
            neighbour_list.update(atoms, cutoff);

            auto [e_pot_l, stress_zz_l] = ducastelle_stress(atoms, neighbour_list, domain.nb_local(), cutoff);

            verlet_step2(atoms.velocities, atoms.forces, dt, m);
            berendsen_thermostat_decomposed(atoms, temp, dt, TAU,domain.nb_local(), m);


            if (i >= Scale) {
                force_local += stress_zz_l(2, 2) / volume(2);
            }

            if (i >= Scale && (i - Scale) % plotEvery == 0) {

                double force_zz_g=MPI::allreduce(force_local / plotEvery,
                                                 MPI_SUM, domain.communicator());
                force_local = 0.0;
                double e_pot_global= MPI::allreduce(e_pot_l, MPI_SUM, domain.communicator());

                domain.disable(atoms);
                remove_mean_velocity(atoms.velocities);

                if (domain.rank() == 0) {
                    std::cout << std::fixed << std::setprecision(15)
                              << temperature_cur(atoms) << "," << force_zz_g << "," << volume(2)
                              << "," << e_pot_global << "," << t << "," << v_z_0 << ","
                              << strainRate << "," << volume(0) << "," << volume(1)
                              << std::endl;
                    results<< std::fixed << std::setprecision(15)
                            << temperature_cur(atoms) << "," << force_zz_g << "," << volume(2)
                            << "," << e_pot_global << "," << t << "," << v_z_0 << ","
                            << strainRate << "," << volume(0) << "," << volume(1)
                            << std::endl;

                    write_xyz(traj, atoms);

                }
                domain.enable(atoms);
                domain.update_ghosts(atoms, cutoff * 2.);
                neighbour_list.update(atoms, cutoff);
                ducastelle(atoms, neighbour_list);
            }

            if (i >= Scale) {
                t += dt;
            }
        }

        if (domain.rank() == 0) {
        }
        }
    }
    MPI_Finalize();
    return 0;
}
