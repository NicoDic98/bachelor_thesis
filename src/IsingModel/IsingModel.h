/**
 * @file       IsingModel.h
 * @brief      Declaration of the Ising Model
 * @author     nico
 * @version    0.0.2
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
     * @brief Load model from \p root
     * @param root Group which holds the model attributes
     */
    [[maybe_unused]] explicit IsingModel(HighFive::Group &root);


    /**
     * @brief Calculates the action for the given field phi
     * @param phi Field
     * @return S(phi) (action)
     */
    double get_action(const VectorX &phi) override;

    /**
     * @brief Sets the value of beta to new_beta
     * @param new_beta Inverse temperature
     */
    void set_beta(double new_beta) override {
        BaseModel::set_beta(new_beta);
        sqrt_beta = sqrt(get_beta());
    }

    /**
     * @brief Calculates the magnetization for the given field phi
     * @param phi Field
     * @return m(phi) (magnetization)
     */
    double get_magnetization(const VectorX &phi) override;
    double get_field_squared(const VectorX &phi) override;


    /**
     * @brief Calculates the magnetization squared for the given field \p phi
     * @param phi Field
     * @return m²(phi) (magnetization squared)
     */
    double get_magnetization_squared(const VectorX &phi) override;

    /**
     * @brief Calculates the average energy per site for the given field \p phi
     * @param phi Field
     * @return E(phi) (average energy per site)
     */
    double get_energy(const VectorX &phi) override;

    /**
     * @brief Calculates the average energy per site squared for the given field \p phi
     * @param phi Field
     * @return e²(phi) (average energy per site squared)
     */
    double get_energy_squared(const VectorX &phi) override;

    /**
     * @brief Checks the dimensions of internal vectors and matrices with regard to the given field phi
     * @param phi Field
     * @return True, if all dimension checks are passed. False, if any dimension check fails.
     */
    [[nodiscard]] bool check_dimensions(const VectorX &phi) const override;

    void ergodicity_jump(VectorX &phi) override;

    /**
     * @brief Returns the coarsent model with respect to the given \p InterpolationType_
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
     * @brief Return the model at \p root
     * @param root Group with model parameters as attributes
     * @return Model loaded with parameters at \p root
     */
    IsingModel *get_model_at(HighFive::Group &root) override;

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

    /**
     * @brief Updates the internal attributes from the \a RootModel
     */
    void pull_attributes_from_root() override;

    /**
     * @brief Returns an empty field, useful for the starting of a Multi Level run
     * @return Empty d.o.f. field
     */
    VectorX get_empty_field() override;

    /**
     * @brief Dump all parameters of the model as attributes to \p root
     * @param root Group to dump to
     */
    void dumpToH5(HighFive::Group &root) override;

    /**
     * @brief Load ensemble from \p root into \p target
     * @param target Target to load to
     * @param root Dataset to load from
     */
    void load_ensemble(std::vector<VectorX> &target, HighFive::DataSet &root) override;

    /**
     * @brief Dump ensemble \p target as dataset to \p root, with the name \p sub_name
     * @param target Ensemble to be dumped
     * @param root Destination group
     * @param sub_name Name to be used for the dataset
     * @return Dataset to which the ensemble was dumped
     */
    HighFive::DataSet dump_ensemble(std::vector<VectorX> &target, HighFive::Group &root, std::string sub_name) override;

    /**
     * @brief Prints the name of the model
     */
    void print_name() override;

    /**
     * @brief Prints the dimensions of the stored vectors and matrices
     */
    [[maybe_unused]] void print_dimensions();

    /**
     * @brief Prints out the connectivity matrix
     */
    void print_connectivity_matrix();

    /**
     * @brief Prints out the interpolation matrix
     */
    [[maybe_unused]] void print_interpolation_matrix();

protected:

    /**
     * @brief Calculates the artificial force for the given field phi
     * @param phi Field
     * @return Artificial force according to artificial hamiltonian
     */
    VectorX get_force(const VectorX &phi) override;

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
     * @brief String to be used as key for \a dimension in H5 files
     */
    static const char *dimension_name;

    /**
     * @brief Neighbour extent
     */
    int neighbour_extent;
    /**
     * @brief String to be used as key for \a neighbour_extent in H5 files
     */
    static const char *neighbour_extent_name;

    /**
     * @brief External field
     */
    VectorX h;
    /**
     * @brief String to be used as key for \a h in H5 files
     */
    static const char *h_name;

    /**
     * @brief Generalization field
     */
    VectorX eta;
    /**
     * @brief String to be used as key for \a eta in H5 files
     */
    static const char *eta_name;

    /**
     * @brief Offset of the connectivity matrix
     */
    double connectivity_offset;
    /**
     * @brief String to be used as key for \a connectivity_offset in H5 files
     */
    static const char *connectivity_offset_name;

    /**
     * @brief Symmetric connectivity matrix
     */
    MatrixX k_sym;
    /**
     * @brief Inverse of \a k_sym
     */
    MatrixX k_sym_inverse;
    /**
     * @brief String to be used as key for \a k_sym in H5 files
     */
    static const char *k_sym_name;
    /**
     * @brief Setter for \a k_sym
     * @param k_sym_new New value for \a k_sym
     */
    void set_k_sym(const MatrixX &k_sym_new);

    /**
     * @brief Asymmetric connectivity matrix
     */
    MatrixX k_rec;
    /**
     * @brief String to be used as key for \a k_rec in H5 files
     */
    static const char *k_rec_name;

    /**
     * @brief Interpolation matrix
     */
    MatrixX InterpolationMatrix;
    /**
     * @brief String to be used as key for \a InterpolationMatrix in H5 files
     */
    static const char *InterpolationMatrix_name;

    /**
     * @brief Reference to the next finer Level in Multi Level mode or the model from which this model is a copy,
     *        otherwise reference to \c *this.
     */
    const IsingModel &RootModel;

    /**
     * @brief Default name of the IsingModel
     */
    static const char *IsingModel_name;

    /**
     * @brief Add identity*\p offset to the connectivity matrix \a k_sym
     * @param offset Offset to be used
     */
    void add_offset_to_connectivity_matrix(double offset);

    /**
     * @brief Checks the dimensions of internal vectors and matrices
     */
    [[nodiscard]] bool check_internal_dimensions() const { return check_dimensions(h); }
};

#endif //BACHELOR_THESIS_ISINGMODEL_H
