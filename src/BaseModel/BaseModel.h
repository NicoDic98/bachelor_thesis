/**
 * @file       BaseModel.h
 * @brief      Base physical models
 * @author     nico
 * @version    0.0.1
 * @date       27.03.22
 */


#ifndef BACHELOR_THESIS_BASEMODEL_H
#define BACHELOR_THESIS_BASEMODEL_H

#include <iostream>
#include <highfive/H5Easy.hpp>
#include <utility>

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
     */
    explicit BaseModel(double beta_, std::string name_ = BaseModel_name) : beta{beta_}, name{std::move(name_)} {}

    /**
     * @brief Copy constructor of BaseModel
     * @param NewModel Model to copy
     */
    BaseModel(const BaseModel<configuration_type> &NewModel) : beta{NewModel.beta}, name{NewModel.name} {}

    explicit BaseModel(HighFive::Group &root, const std::string &default_name_ = BaseModel_name)
             :beta{}{
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
    virtual inline double get_beta() const;

    /**
     * @brief Sets the value of beta to new_beta
     * @param new_beta Inverse temperature
     */
    virtual inline void set_beta(double new_beta);

    /**
     * @brief Calculates the magnetization for the given field phi
     * @param phi Field
     * @return m(phi) (magnetization)
     */
    virtual double get_magnetization(const configuration_type &phi) { return 0; }


    /**
     * @brief Checks the dimensions of internal vectors and matrices with regard to the given field phi
     * @param phi Field
     * @return True, if all dimension checks are passed. False, if any dimension check fails.
     */
    virtual bool check_dimensions(const configuration_type &phi) const = 0;

    /**
     * @brief Updates the momentum pi and returns the new momentum pi_new
     * @param phi
     * @param pi
     * @param step_size
     * @return pi_new
     */
    virtual void update_pi(configuration_type &phi, configuration_type &pi, double step_size);

    /**
     * @brief Updates the field phi and returns the new field phi_new
     * @param phi
     * @param pi
     * @param step_size
     * @return phi_new
     */
    virtual void update_phi(configuration_type &phi, configuration_type &pi, double step_size);

    /**
     * @brief Returns a configuration using the standard constructor of the configuration type
     * @return Sample object, which represents a configuration
     */
    virtual configuration_type get_dof_sample();

    /**
     * @brief Returns the coarsent model with respect to the given interpolation matrix
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

    virtual void pull_attributes_from_finer_level() = 0;

    /**
     * @brief Returns an empty field, useful for the starting of a Multi Level run
     * @return Empty field d.o.f. field
     */
    virtual configuration_type get_empty_field() = 0;

    virtual void dumpToH5(HighFive::Group &root);

    virtual void
    load_ensemble(std::vector<configuration_type> &target, HighFive::DataSet &root) = 0;

    /**
     * @brief Prints the name of the model
     */
    virtual void print_name();

protected:

    /**
     * @brief Calculates the artificial force for the given field phi
     * @param phi Field
     * @return Artificial force according to artificial hamiltonian
     */
    virtual configuration_type get_force(const configuration_type &phi) = 0;

    static const char *beta_name;
    static const char *model_name_key;
    static const char *BaseModel_name;

    std::string name;

private:
    /**
     * @brief Inverse temperature
     */
    double beta;
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
void BaseModel<configuration_type>::update_pi(configuration_type &phi, configuration_type &pi, double step_size) {
    pi += step_size * get_force(phi);
}

template<class configuration_type>
void BaseModel<configuration_type>::update_phi(configuration_type &phi, configuration_type &pi, double step_size) {
    phi += step_size * pi;
}

template<class configuration_type>
configuration_type BaseModel<configuration_type>::get_dof_sample() {
    return configuration_type();
}

template<class configuration_type>
void BaseModel<configuration_type>::print_name() {
    std::cout << name << std::endl;
}

template<class configuration_type>
void BaseModel<configuration_type>::dumpToH5(HighFive::Group &root) {
    if (root.hasAttribute(beta_name)) {
        HighFive::Attribute temp = root.getAttribute(beta_name);
        temp.write(beta);//type should always be double, no size checks
    } else {
        root.createAttribute(beta_name, beta);
    }
    if (root.hasAttribute(model_name_key)) {
        HighFive::Attribute temp = root.getAttribute(model_name_key);
        temp.write(name);
    } else {
        root.createAttribute(model_name_key, name);
    }
}


#endif //BACHELOR_THESIS_BASEMODEL_H
