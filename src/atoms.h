//
// Created by walid on 7/3/24.
//

#ifndef YAMD_ATOMS_H
#define YAMD_ATOMS_H

#include "types.h"



class Atoms {
  public:
    Names_t names;
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Masses_t masses;



    Atoms(const Positions_t &p) :
          positions{p}, velocities{3, p.cols()}, forces{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }

    Atoms(const Positions_t &p, const Velocities_t &v) :
          positions{p}, velocities{v}, forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.setOnes();
    }



    Atoms(const Names_t& names, const Positions_t& positions)
        : names(names), positions(positions), velocities(3, positions.cols()), forces(3, positions.cols()), masses(positions.cols()) {
        velocities.setZero();
        forces.setZero();
        for (auto i{0}; i < positions.cols(); i++) {
            masses(i) = ELEM_NAME_TO_MASS.at(names[i]);
        }
    }



    Atoms(int nb_atoms)
        // Initialize the object with nb_atoms
        : positions{3, nb_atoms},
          velocities{3, nb_atoms},
          forces{3, nb_atoms}
          {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }

    void resize(size_t new_size) {
        positions.conservativeResize(3, new_size);
        velocities.conservativeResize(3, new_size);
        forces.conservativeResize(3, new_size);
        masses.conservativeResize( new_size);

    }

    double kinetic_energy() const {
        // use `colwise()` and `transpose()` to perform the sum using only Eigen
        // functions
        return ((Eigen::ArrayXd)(velocities.colwise().squaredNorm().transpose() *
                                 masses * 0.5))
            .sum();
    }

    double kinetic_energy(int nb_local, double mass) const {
        double total_ke = 0.0;  // Initialize total kinetic energy

        // Loop over each atom up to nb_local
        for (int i = 0; i < nb_local; ++i) {
            double vec_norm = 0.0;  // Initialize the norm of the velocity vector

            // Compute the squared norm of the i-th velocity vector (3 components: x, y, z)
            for (int j = 0; j < velocities.rows(); ++j) {
                vec_norm += velocities(j, i) * velocities(j, i);
            }

            // Multiply the squared norm by mass/2 and add to the total kinetic energy
            total_ke += 0.5 * mass * vec_norm;
        }

        return total_ke;  // Return the total kinetic energy
    }


};


#endif // YAMD_ATOMS_H
