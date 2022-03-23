/**
 * @file       IsingModel.h
 * @brief      Declaration of the Ising Model Including its Connectivity matrix and Force term in the artificial hamiltonian
 * @author     nico
 * @version    0.0.1
 * @date       23.03.22
 */


#ifndef BACHELOR_THESIS_ISINGMODEL_H
#define BACHELOR_THESIS_ISINGMODEL_H

#include <eigen3/Eigen/Dense>
#include <vector>

namespace {
    /**
     * @brief Vector type definition
     */
    typedef Eigen::VectorXd VectorX;
    typedef Eigen::MatrixXd MatrixX;
}

class IsingModel {
public:
    IsingModel(double beta_, VectorX h_, VectorX eta_, double offset_, int dimension_, int neighbour_extent_,
               int grid_size_);

    IsingModel(double beta_, VectorX h_, VectorX eta_, double offset_, MatrixX k_sym_, MatrixX k_rec_);

    double get_action(const VectorX &phi);

    VectorX get_force(const VectorX &phi);

    double get_magnetization(const VectorX &phi);

    void print_connectivity_matrix();

private:
    double beta;
    double sqrt_beta;
    VectorX h;
    VectorX eta;
    double offset;
    MatrixX k_sym;
    MatrixX k_rec;

    void fill_connectivity_matrix(int dimension, int neighbour_extent, int grid_size);

    void add_offset_to_connectivity_matrix();

};

inline int int_pow(int x, int y) {
    assert(y > 0);
    int result = 1;
    for (int i = 0; i < y; ++i) {
        result *= x;
    }
    return result;
}

#endif //BACHELOR_THESIS_ISINGMODEL_H
