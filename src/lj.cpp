#include "lj.h"
#include <cmath>
#include "lj_direct_summation.h"


// Function to compute the Lennard-Jones potential
double lj_cutoff_potential(double r, double epsilon, double sigma, double cutoff) {
    if (r > cutoff) {
        return 0.0;
    }
    return lj_potential(r,epsilon,sigma) - lj_potential(cutoff,epsilon,sigma);
}


double lj_cutoff_direct_summation(Atoms &atoms, double epsilon, double sigma, double cutoff) {
    double potential_energy = 0.0;
    atoms.forces.setZero();
    size_t nb_atoms = atoms.nb_atoms();

    NeighborList neighborList;
    auto [seed, neighbors] = neighborList.update(atoms, cutoff);

    for (size_t i = 0; i < nb_atoms; ++i) {
        for (size_t j = seed[i]; j < seed[i + 1]; ++j) {
            size_t neighbor_index = neighbors[j];
            if (i < neighbor_index) {
                Eigen::Vector3d r_ij =
                    atoms.positions.col(i).matrix() -
                    atoms.positions.col(neighbor_index).matrix();
                double r = r_ij.norm();
                if (r > cutoff) continue;

                Eigen::Vector3d r_hat_ij = r_ij / r;

                double potential = lj_cutoff_potential(r, epsilon, sigma, cutoff);
                double potential_derivative =
                    lj_potential_derivative(r, epsilon, sigma);

                potential_energy += potential;

                Eigen::Vector3d force = r_hat_ij * potential_derivative;
                atoms.forces.col(i) += force.array();
                atoms.forces.col(neighbor_index) -= force.array();
            }
        }
    }
    return potential_energy;
}