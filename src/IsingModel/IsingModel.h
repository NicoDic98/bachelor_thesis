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
#include <BaseModel.h>

/**
 * @brief Definition of the Ising model.
 */
class IsingModel : public BaseModel<VectorX> {
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
     * @brief Coarsening constructor of IsingModel
     * @param NewModel Finer model
     * @param InterpolationMatrix Interpolation matrix to use for the coarsening
     */
    IsingModel(const IsingModel &NewModel, const MatrixX &InterpolationMatrix);

    /**
     * @brief Copy constructor of IsingModel
     * @param NewModel Model to copy
     */
    IsingModel(const IsingModel &NewModel);


    /**
     * @brief Calculates the action for the given field phi
     * @param phi Field
     * @return S(phi) (action)
     */
    double get_action(const VectorX &phi) override;

    /**
     * @brief Calculates the magnetization for the given field phi
     * @param phi Field
     * @return m(phi) (magnetization)
     */
    double get_magnetization(const VectorX &phi) override;

    /**
     * @brief Prints out the connectivity matrix
     */
    void print_connectivity_matrix();

    /**
     * @brief Sets the value of beta to new_beta
     * @param new_beta Inverse temperature
     */
    void set_beta(double new_beta) override {
        BaseModel::set_beta(new_beta);
        sqrt_beta = sqrt(get_beta());
    }

    /**
     * @brief Checks the dimensions of internal vectors and matrices with regard to the given field phi
     * @param phi Field
     * @return True, if all dimension checks are passed. False, if any dimension check fails.
     */
    bool check_dimensions(const VectorX &phi) override;

    /**
     * @brief Prints the name of the model
     */
    void print_name() override;

    /**
     * @brief Returns the coarsent model with respect to the given interpolation matrix
     * @param InterpolationMatrix Interpolation matrix to use for the coarsening
     * @return Coursed model
     */
    IsingModel *get_coarser_model(const MatrixX &InterpolationMatrix) override;

    /**
     * @brief Return a copy of the model
     * @return Copy of the model
     */
    IsingModel *get_copy_of_model() override;

private:
    /**
     * @brief sqrt(Inverse temperature)
     */
    double sqrt_beta;

    /**
     * @brief Dimension of the hypercube
     */
    int dimension;

    /**
     * @brief External field
     */
    VectorX h;

    /**
     * @brief Generalization field
     */
    VectorX eta;

    /**
     * @brief Symmetric connectivity matrix
     */
    MatrixX k_sym;

    /**
     * @brief Asymmetric connectivity matrix
     */
    MatrixX k_rec;

    /**
     * @brief Fills the connectivity matrix \a k_sym for the given hyper cube of dimension \a dimension and side length \p grid_size
     * @param neighbour_extent Range at which nearest neighbours interact (along the axis)
     * @param grid_size Side length
     */
    void fill_connectivity_matrix(int neighbour_extent, int grid_size);

    /**
     * @brief Add identity*\p offset to the connectivity matrix \a k_sym
     * @param offset Offset to be used
     */
    void add_offset_to_connectivity_matrix(double offset);

    /**
     * @brief Checks the dimensions of internal vectors and matrices
     */
    bool check_internal_dimensions() { return check_dimensions(h); }

protected:

    /**
     * @brief Calculates the artificial force for the given field phi
     * @param phi Field
     * @return Artificial force according to artificial hamiltonian
     */
    VectorX get_force(const VectorX &phi) override;

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
