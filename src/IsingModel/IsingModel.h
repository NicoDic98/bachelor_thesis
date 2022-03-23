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
    explicit IsingModel(int dimension_, int grid_size_, int neighbour_extent_, double offset);

    double get_action(const VectorX &phi);

    VectorX get_force(const VectorX &phi);

    double get_magnetization(const VectorX &phi);

    void print_connectivity_matrix();

private:
    int dimension;
    int grid_size;
    int neighbour_extent;
    MatrixX connectivity_matrix;

    void fill_connectivity_matrix();
    void add_offset_to_connectivity_matrix(double offset);
    void generate_index_offsets(std::vector<long> &index_offsets, long size);
};


#endif //BACHELOR_THESIS_ISINGMODEL_H
