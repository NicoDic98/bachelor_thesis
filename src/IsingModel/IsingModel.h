/**
 * @file       IsingModel.h
 * @brief      Declaration of the Ising Model Including its Connectivity matrix and Force term in the artificial hamiltonian
 * @author     nico
 * @version    0.0.1
 * @date       23.03.22
 * @todo       Use sparse matrices for the connectivity matrices
 */


#ifndef BACHELOR_THESIS_ISINGMODEL_H
#define BACHELOR_THESIS_ISINGMODEL_H


#include <vector>
#include <MyTypes.h>

/**
 * @brief Definition of the Ising model.
 */
class IsingModel {
public:
    /**
     * @brief Constructor of IsingModel, with generation of connectivity matrix
     * @param beta_ Inverse temperature
     * @param h_ External field
     * @param eta_ Generalisation field
     * @param offset_ Offset used to shift the connectivity matrix by identity*offset_
     * @param dimension_ Dimension of the Ising model
     * @param neighbour_extent_ Range of coupling spins
     * @param grid_size_ Spatial extend of the square/cubic/... grid
     */
    IsingModel(double beta_, VectorX h_, VectorX eta_, double offset_, int dimension_, int neighbour_extent_,
               int grid_size_);

    /**
     * @brief Constructor of IsingModel, expecting already calculated k_sym and k_rec
     * @param beta_ Inverse temperature
     * @param h_ External field
     * @param eta_ Generalisation field
     * @param offset_ Offset used to shift the connectivity matrix by identity*offset_
     * @param k_sym_ Square connectivity matrix
     * @param k_rec_ Rectangular connectivity matrix
     */
    IsingModel(double beta_, VectorX h_, VectorX eta_, double offset_, MatrixX k_sym_, MatrixX k_rec_);

    /**
     * @brief Calculates the action for the given field phi
     * @param phi Field
     * @return S(phi) (action)
     */
    double get_action(const VectorX &phi);

    /**
     * @brief Calculates the artificial force for the given field phi
     * @param phi Field
     * @return Artificial force according to artificial hamiltonian
     */
    VectorX get_force(const VectorX &phi);

    /**
     * @brief Calculates the magnetization for the given field phi
     * @param phi Field
     * @return m(phi) (magnetization)
     */
    double get_magnetization(const VectorX &phi);

    /**
     * @brief Prints out the connectivity matrix
     */
    void print_connectivity_matrix();

    /**
     * @brief Sets the value of beta to new_beta
     * @param new_beta Inverse temperature
     */
    void set_beta(double new_beta) {
        beta = new_beta;
        sqrt_beta = sqrt(beta);
    }

    /**
     * @brief Returns beta
     * @return Inverse temperature
     */
    double get_beta() const { return beta; }

    bool check_dimensions(const VectorX &phi);

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

    bool check_internal_dimensions() { return check_dimensions(h); }
};

/**
 * @brief Simple Power function to calculate the powers of integers
 * @param x Base
 * @param y Exponent
 * @return Base^Exponent
 */
inline int int_pow(int x, int y) {
    assert(y > 0);
    int result = 1;
    for (int i = 0; i < y; ++i) {
        result *= x;
    }
    return result;
}

#endif //BACHELOR_THESIS_ISINGMODEL_H
