/**
 * @file       XYModel.h
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       18.05.22
 */


#ifndef MULTILEVELHMC_XYMODEL_H
#define MULTILEVELHMC_XYMODEL_H

#include <BaseModel.h>
#include <vector>
#include <MyTypes.h>

/**
 * @brief Definition of the XY model
 */
class XYModel : public BaseModel<MultiVectorX> {
public:
    XYModel(double beta_,double eta_, MultiVectorX h_, int dimension_, int neighbour_extent_, int grid_size_);

    /**
     * @brief Coarsening constructor of XYModel
     * @param NewModel Finer model
     * @param InterpolationType_ Interpolation type to use for the coarsening
     */
    XYModel(const XYModel &NewModel, InterpolationType InterpolationType_);

    /**
     * @brief Copy constructor of XYModel
     * @param NewModel Model to copy
     */
    XYModel(const XYModel &NewModel);

    /**
     * @brief Load model from \p root
     * @param root Group which holds the model attributes
     */
    [[maybe_unused]] explicit XYModel(HighFive::Group &root);

    /**
     * @brief Calculates the action for the given field phi
     * @param phi Field
     * @return S(phi) (action)
     */
    double get_action(const MultiVectorX &phi) override;

    /**
     * @brief Calculates the artificial energy for the given field phi
     * @param phi Field
     * @return H(phi) (artificial energy)
     */
    double get_artificial_energy(const MultiVectorX &phi, const MultiVectorX &pi) override;

    MultiVectorX get_pi(std::default_random_engine &generator) override;

    double get_energy(const MultiVectorX &phi) override;

    double get_vector_length_squared(const MultiVectorX &phi)override;

    /**
     * @brief Checks the dimensions of internal vectors and matrices with regard to the given field phi
     * @param phi Field
     * @return True, if all dimension checks are passed. False, if any dimension check fails.
     */
    [[nodiscard]] bool check_dimensions(const MultiVectorX &phi) const override;

    /**
     * @brief Updates the momentum \p pi in place
     * @param phi
     * @param pi
     * @param step_size
     */
    void update_pi(MultiVectorX &phi, MultiVectorX &pi, double step_size) override;

    /**
     * @brief Updates the field \p phi in place
     * @param phi
     * @param pi
     * @param step_size
     */
    void update_phi(MultiVectorX &phi, MultiVectorX &pi, double step_size) override;

    /**
     * @brief Returns the coarsent model with respect to the given \p InterpolationType_
     * @param InterpolationType_ Interpolation type to use for the coarsening
     * @return Coursed model
     */
    XYModel *get_coarser_model(InterpolationType InterpolationType_) override;

    /**
     * @brief Return a copy of the model
     * @return Copy of the model
     */
    XYModel *get_copy_of_model() override;


    /**
     * @brief Return the model at \p root
     * @param root Group with model parameters as attributes
     * @return Model loaded with parameters at \p root
     */
    XYModel *get_model_at(HighFive::Group &root) override;

    /**
     * @brief Update internal Fields with the d.o.f. field \p phi of the finer level.
     * @param phi d.o.f. field
     */
    void update_fields(const MultiVectorX &phi) override;

    /**
     * @brief Updates finer field \p phia using the \a InterpolationMatrix and the coarser field \p phi2a
     * @param phi2a Coarse field
     * @param phia Fine field
     */
    void interpolate(const MultiVectorX &phi2a, MultiVectorX &phia) override;

    /**
     * @brief Updates the internal attributes from the \a RootModel
     */
    void pull_attributes_from_root() override;

    /**
     * @brief Returns an empty field, useful for the starting of a Multi Level run
     * @return Empty d.o.f. field
     */
    MultiVectorX get_empty_field() override;

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
    void load_ensemble(std::vector<MultiVectorX> &target, HighFive::DataSet &root) override;

    /**
     * @brief Dump ensemble \p target as dataset to \p root, with the name \p sub_name
     * @param target Ensemble to be dumped
     * @param root Destination group
     * @param sub_name Name to be used for the dataset
     * @return Dataset to which the ensemble was dumped
     */
    HighFive::DataSet dump_ensemble(std::vector<MultiVectorX> &target,
                                    HighFive::Group &root, std::string sub_name) override;

    void renormalize(MultiVectorX &phi) override;


protected:

    /**
     * @brief Calculates the artificial force for the given field phi
     * @param phi Field
     * @return Artificial force according to artificial hamiltonian
     */
    MultiVectorX get_force(const MultiVectorX &phi) override;

private:
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

    double eta;

    static const char *eta_name;

    /**
     * @brief External field
     */
    MultiVectorX h;
    /**
     * @brief String to be used as key for \a h in H5 files
     */
    static const char *h_name;

    /**
     * @brief Symmetric connectivity matrix
     */
    MatrixX k_sym;
    /**
     * @brief String to be used as key for \a k_sym in H5 files
     */
    static const char *k_sym_name;

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
    const XYModel &RootModel;

    /**
     * @brief Default name of the XYModel
     */
    static const char *XYModel_name;

    /**
     * @brief Checks the dimensions of internal vectors and matrices
     */
    [[nodiscard]] bool check_internal_dimensions() const { return check_dimensions(h); }
};


#endif //MULTILEVELHMC_XYMODEL_H
