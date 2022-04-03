/**
 * @file       IsingModel.h
 * @brief      Declaration of the Ising Model
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
     * @param offset_ Offset used to shift the connectivity matrix by \c identity*(offset_+ needed_offset)
     * @param dimension_ Dimension of the Ising model
     * @param neighbour_extent_ Range of coupling spins
     * @param grid_size_ Spatial extend of the square/cubic/... grid
     */
    IsingModel(double beta_, VectorX h_, VectorX eta_, double offset_, int dimension_, int neighbour_extent_,
               int grid_size_);

    /**
     * @brief Coarsening constructor of IsingModel
     * @param NewModel Finer model
     * @param InterpolationType_ Interpolation type to use for the coarsening
     */
    IsingModel(const IsingModel &NewModel, InterpolationType InterpolationType_);

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
     * @brief Prints out the interpolation matrix
     */
    [[maybe_unused]] void print_interpolation_matrix();

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
    [[nodiscard]] bool check_dimensions(const VectorX &phi) const override;

    /**
     * @brief Prints the name of the model
     */
    void print_name() override;

    /**
     * @brief Prints the dimensions of the stored vectors and matrices
     */
    [[maybe_unused]] void print_dimensions();

    /**
     * @brief Returns the coarsent model with respect to the given interpolation matrix
     * @param InterpolationType_ Interpolation type to use for the coarsening
     * @return Coursed model
     */
    IsingModel *get_coarser_model(InterpolationType InterpolationType_) override;

    /**
     * @brief Return a copy of the model
     * @return Copy of the model
     */
    IsingModel *get_copy_of_model() override;

    /**
     * @brief Update internal Fields with the d.o.f. field \p phi of the finer level.
     * @param phi d.o.f. field
     */
    void update_fields(const VectorX &phi) override;

    /**
     * @brief Updates finer field \p phia using the \a InterpolationMatrix and the coarser field \p phi2a
     * @param phi2a Coarse field
     * @param phia Fine field
     */
    void interpolate(const VectorX &phi2a, VectorX &phia) override;

    void pull_attributes_from_finer_level() override;

    /**
     * @brief Returns an empty field, useful for the starting of a Multi Level run
     * @return Empty field d.o.f. field
     */
    VectorX get_empty_field() override;

    /**
     * @brief Dumps all data as attributes to the H5 \p file at \p path
     * @param file File to dump to
     * @param path Path to dump to
     */
    void dumpToH5(HighFive::File &file, std::string path) override;

private:
    /**
     * @brief sqrt(Inverse temperature)
     */
    double sqrt_beta;

    /**
     * @brief Dimension of the hypercube
     */
    int dimension;
    static const char *dimension_name;

    /**
     * @brief Side length of the hypercube
     */
    int grid_side_length;
    static const char *grid_side_length_name;

    /**
     * @brief External field
     */
    VectorX h;
    static const char *h_name;

    /**
     * @brief Generalization field
     */
    VectorX eta;
    static const char *eta_name;

    /**
     * @brief Symmetric connectivity matrix
     */
    MatrixX k_sym;
    MatrixX k_sym_inverse;
    static const char *k_sym_name;

    /**
     * @brief Asymmetric connectivity matrix
     */
    MatrixX k_rec;
    static const char *k_rec_name;

    /**
     * @brief Interpolation matrix
     */
    MatrixX InterpolationMatrix;
    static const char *InterpolationMatrix_name;

    /**
     * @brief Reference to the next finer Level in Multi Level mode, otherwise reference to \c *this.
     */
    const IsingModel &FinerModel;

    /**
     * @brief Fills the connectivity matrix \a k_sym for the given hyper cube of dimension \a dimension and side length \p grid_size
     * @param neighbour_extent Range at which nearest neighbours interact (along the axis)
     * @param grid_size Side length
     */
    void fill_connectivity_matrix(int neighbour_extent, int grid_size);

    /**
     * @brief Fills the \a InterpolationMatrix
     * @param InterpolationType_ Type of interpolation to be used
     * @param fine_size Finer grid total size
     * @param fine_grid_side_length Finer grid hyper cube side length
     * @return Coarse grid side length
     */
    int fill_interpolation_matrix(InterpolationType InterpolationType_, long fine_size, int fine_grid_side_length);

    /**
     * @brief Add identity*\p offset to the connectivity matrix \a k_sym
     * @param offset Offset to be used
     */
    void add_offset_to_connectivity_matrix(double offset);

    /**
     * @brief Checks the dimensions of internal vectors and matrices
     */
    [[nodiscard]] bool check_internal_dimensions() const { return check_dimensions(h); }

    void set_k_sym(const MatrixX &k_sym_new);

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
    assert(y >= 0);
    int result = 1;
    for (int i = 0; i < y; ++i) {
        result *= x;
    }
    return result;
}

#endif //BACHELOR_THESIS_ISINGMODEL_H
