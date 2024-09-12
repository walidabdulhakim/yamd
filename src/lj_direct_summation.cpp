#include "lj_direct_summation.h"
#include <cmath>


// Function to compute the Lennard-Jones potential
double lj_potential(double r, double epsilon, double sigma) {
    double r6 = std::pow(sigma / r, 6);
    double r12 = r6 * r6;
    return 4 * epsilon * (r12 - r6);
}

// Function to compute the derivative of the Lennard-Jones potential with respect to distance
double lj_potential_derivative(double r, double epsilon, double sigma) {
    double r6 = std::pow(sigma / r, 6);
    double r12 = r6 * r6;
    return 24 * epsilon * (2 * r12 - r6) / r;
}

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    double potential_energy = 0.0;
    atoms.forces.setZero();
    size_t nb_atoms = atoms.nb_atoms();

    for (size_t i = 0; i < nb_atoms; ++i) {
        for (size_t j = i + 1; j < nb_atoms; ++j) {
            Eigen::Vector3d r_ij = atoms.positions.col(i).matrix() - atoms.positions.col(j).matrix();
            double r = r_ij.norm();
            Eigen::Vector3d r_hat_ij = r_ij / r;

            // Use the modular potential and derivative functions
            double potential = lj_potential(r, epsilon, sigma);
            double potential_derivative = lj_potential_derivative(r, epsilon, sigma);

            potential_energy += potential;

            Eigen::Vector3d force = r_hat_ij * potential_derivative;
            atoms.forces.col(i) += force.array();
            atoms.forces.col(j) -= force.array();
        }
    }
    return potential_energy;
}

