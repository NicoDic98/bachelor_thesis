/**
 * @file       BaseModel.h
 * @brief      Base physical models
 * @author     nico
 * @version    0.0.2
 * @date       27.03.22
 */


#ifndef BACHELOR_THESIS_BASEMODEL_H
#define BACHELOR_THESIS_BASEMODEL_H

#include <iostream>
#include <highfive/H5File.hpp>
#include <utility>
#include <MyTypes.h>

/**
 * @brief Template for the Abstract Class BaseModel
 * @tparam configuration_type Datatype, which is used for the configurations of the model
 */
template<class configuration_type>
class BaseModel {
public:

    /**
     * @brief Most Basic constructor of BaseModel
     * @param beta_ Inverse temperature of the model
     * @param name_ Name of the Model, default="BaseModel"
     */
    explicit BaseModel(double beta_, std::string name_ = BaseModel_name) : beta{beta_}, name{std::move(name_)} {}

    /**
     * @brief Copy constructor of BaseModel
     * @param NewModel Model to copy
     */
    BaseModel(const BaseModel<configuration_type> &NewModel) : beta{NewModel.beta}, name{NewModel.name} {}

    /**
     * @brief Load model from \p root
     * @param root Group which holds the model attributes
     * @param default_name_ Default name to be used, when there is no model name present in \p root
     */
    explicit BaseModel(HighFive::Group &root, const std::string &default_name_ = BaseModel_name)
            : beta{} {
        root.getAttribute(beta_name).read(beta);
        if (root.hasAttribute(model_name_key)) {
            root.getAttribute(model_name_key).read(name);
        } else {
            name = default_name_;
        }
    }

    /**
     * @brief Destructor of BaseModel
     */
    virtual ~BaseModel() = default;

    /**
     * @brief Calculates the action for the given field phi
     * @param phi Field
     * @return S(phi) (action)
     */
    virtual double get_action(const configuration_type &phi) = 0;

    /**
     * @brief Returns beta
     * @return Inverse temperature
     */
    [[nodiscard]] virtual inline double get_beta() const;

    /**
     * @brief Sets the value of beta to new_beta
     * @param new_beta Inverse temperature
     */
    virtual inline void set_beta(double new_beta);

    /**
     * @brief Calculates the magnetization for the given field \p phi
     * @param phi Field
     * @return m(phi) (magnetization)
     */
    [[maybe_unused]] virtual double get_magnetization(const configuration_type &phi) { return 0; }

    [[maybe_unused]] virtual double get_field_squared(const configuration_type &phi) { return 0; }

    /**
     * @brief Calculates the magnetization squared for the given field \p phi
     * @param phi Field
     * @return m²(phi) (magnetization squared)
     */
    [[maybe_unused]] virtual double get_magnetization_squared(const configuration_type &phi) { return 0; }

    /**
     * @brief Calculates the average energy per site for the given field \p phi
     * @param phi Field
     * @return E(phi) (average energy per site)
     */
    [[maybe_unused]] virtual double get_energy(const configuration_type &phi) { return 0; }

    /**
     * @brief Calculates the average energy per site squared for the given field \p phi
     * @param phi Field
     * @return e²(phi) (average energy per site squared)
     */
    [[maybe_unused]] virtual double get_energy_squared(const configuration_type &phi) { return 0; }


    /**
     * @brief Checks the dimensions of internal vectors and matrices with regard to the given field \p phi
     * @param phi Field
     * @return True, if all dimension checks are passed. False, if any dimension check fails.
     */
    virtual bool check_dimensions(const configuration_type &phi) const = 0;

    /**
     * @brief Updates the momentum \p pi in place
     * @param phi
     * @param pi
     * @param step_size
     */
    virtual void update_pi(configuration_type &phi, configuration_type &pi, double step_size)=0;

    /**
     * @brief Updates the field \p phi in place
     * @param phi
     * @param pi
     * @param step_size
     */
    virtual void update_phi(configuration_type &phi, configuration_type &pi, double step_size)=0;

    virtual void ergodicity_jump(configuration_type &phi);

    /**
     * @brief Returns a configuration using the standard constructor of the configuration type
     * @return Sample object, which represents a configuration
     */
    [[maybe_unused]] [[maybe_unused]] virtual configuration_type get_dof_sample();

    /**
     * @brief Returns the coarsent model with respect to the given \p InterpolationType_
     * @param InterpolationType_ Interpolation type to use for the coarsening
     * @return Coursed model
     */
    virtual BaseModel<configuration_type> *get_coarser_model(InterpolationType InterpolationType_) = 0;

    /**
     * @brief Return a copy of the model
     * @return Copy of the model
     */
    virtual BaseModel<configuration_type> *get_copy_of_model() = 0;

    /**
     * @brief Return the model at \p root
     * @param root Group with model parameters as attributes
     * @return Model loaded with parameters at \p root
     */
    virtual BaseModel<configuration_type> *get_model_at(HighFive::Group &root) = 0;

    /**
     * @brief Update internal Fields with the d.o.f. field \p phi of the finer level
     * @param phi d.o.f. field
     */
    virtual void update_fields(const configuration_type &phi) = 0;

    /**
     * @brief Updates finer field \p phia using the interpolation of the coarser field \p phi2a
     * @param phi2a Coarse field
     * @param phia Fine field
     */
    virtual void interpolate(const configuration_type &phi2a, configuration_type &phia) = 0;

    /**
     * @brief Updates the internal attributes from the finer level model
     */
    virtual void pull_attributes_from_root() = 0;

    /**
     * @brief Returns an empty field, useful for the starting of a Multi Level run
     * @return Empty d.o.f. field
     */
    virtual configuration_type get_empty_field() = 0;

    /**
     * @brief Dump all parameters of the model as attributes to \p root
     * @param root Group to dump to
     */
    virtual void dumpToH5(HighFive::Group &root);

    /**
     * @brief Load ensemble from \p root into \p target
     * @param target Target to load to
     * @param root Dataset to load from
     */
    virtual void
    load_ensemble(std::vector<configuration_type> &target, HighFive::DataSet &root) = 0;

    /**
     * @brief Dump ensemble \p target as dataset to \p root, with the name \p sub_name
     * @param target Ensemble to be dumped
     * @param root Destination group
     * @param sub_name Name to be used for the dataset
     * @return Dataset to which the ensemble was dumped
     */
    virtual HighFive::DataSet
    dump_ensemble(std::vector<configuration_type> &target, HighFive::Group &root, std::string sub_name) = 0;

    /**
     * @brief Prints the name of the model
     */
    [[maybe_unused]] virtual void print_name();

protected:

    /**
     * @brief Calculates the artificial force for the given field \p phi
     * @param phi Field
     * @return Artificial force according to artificial hamiltonian
     */
    virtual configuration_type get_force(const configuration_type &phi) = 0;

    /**
     * @brief Name of the current object
     */
    std::string name;

private:
    /**
     * @brief Inverse temperature
     */
    double beta;
    /**
     * @brief String to be used as key for \a beta in H5 files
     */
    [[maybe_unused]] static const char *beta_name;

    /**
     * @brief String to be used as key for the model \a name in H5 files
     */
    [[maybe_unused]] static const char *model_name_key;

    /**
     * @brief Default name of the BaseModel
     */
    [[maybe_unused]] static const char *BaseModel_name;
};

template<class configuration_type> const char *BaseModel<configuration_type>::beta_name{"beta"};
template<class configuration_type> const char *BaseModel<configuration_type>::model_name_key{"model"};
template<class configuration_type> const char *BaseModel<configuration_type>::BaseModel_name{"BaseModel"};

template<class configuration_type>
inline double BaseModel<configuration_type>::get_beta() const {
    return beta;
}

template<class configuration_type>
inline void BaseModel<configuration_type>::set_beta(double new_beta) {
    beta = new_beta;
}

template<class configuration_type>
[[maybe_unused]] configuration_type BaseModel<configuration_type>::get_dof_sample() {
    return configuration_type();
}

template<class configuration_type>
[[maybe_unused]] void BaseModel<configuration_type>::print_name() {
    std::cout << name << std::endl;
}

template<class configuration_type>
void BaseModel<configuration_type>::dumpToH5(HighFive::Group &root) {
    write_static_size(beta, root, beta_name);
    if (root.hasAttribute(model_name_key)) {
        root.deleteAttribute(model_name_key);
        //TODO maybe there is a nicer fix for this
    }
    root.createAttribute(model_name_key, name);
}

template<class configuration_type>
void BaseModel<configuration_type>::ergodicity_jump(configuration_type &phi) {
    //general ergodicity jump is to do nothing
}


#endif //BACHELOR_THESIS_BASEMODEL_H
