/**
 * @file       HMCGenerator.h
 * @brief      Declarations of standard HMC
 * @author     nico
 * @version    0.0.2
 * @date       24.03.22
 */


#ifndef BACHELOR_THESIS_HMCGENERATOR_H
#define BACHELOR_THESIS_HMCGENERATOR_H

#include <BaseModel.h>
#include <LeapFrogIntegrator.h>
#include <iostream>
#include <random>

/**
 * @brief Template of the standard HMC algorithm
 * @tparam configuration_type Datatype, which is used for the configurations of the model
 */
template<class configuration_type>
class HMCGenerator {
public:
    /**
     * @brief Standard constructor of the HMC generator
     * @param model_ Model for which to generate ensembles
     * @param amount_of_steps_ Amount of steps to be used in the integration process
     * @param step_size_ Step size in the integration process
     * @param generator_ Random number generator to be used for the HMC
     */
    HMCGenerator(BaseModel<configuration_type> &model_, size_t amount_of_steps_, double step_size_,
                 std::default_random_engine &generator_);

    /**
     * @brief Loads the HMCGenerator from \ p root
     * @param model_ Model for which to generate ensembles
     * @param root Group which has the HMCGenerator parameters as Attributes and ensemble as dataset
     * @param generator_ Random number generator to be used for the HMC
     */
    HMCGenerator(BaseModel<configuration_type> &model_, HighFive::Group &root,
                 std::default_random_engine &generator_);

    /**
     * @brief Loads the HMCGenerator from \ p root
     * @param model_ Model for which to generate ensembles
     * @param root Group which has the HMCGenerator parameters as Attributes and ensemble as dataset
     * @param amount_of_steps_ Amount of steps to be used in the integration process
     * @param step_size_ Step size in the integration process
     * @param generator_ Random number generator to be used for the HMC
     */
    HMCGenerator(BaseModel<configuration_type> &model_, HighFive::Group &root,
                 size_t amount_of_steps_, double step_size_,
                 std::default_random_engine &generator_);

    /**
     * @brief Generate \p amount_of_samples amount of ensembles, starting from \p phiStart and doing
     *        \p amount_of_thermalization_steps thermalization steps in advance
     * @param phiStart Starting field
     * @param amount_of_samples Amount of samples to take
     * @param amount_of_thermalization_steps Amount of thermalization steps
     * @param expand Rather to expand the internal \a ensembles vector or not
     * @return Acceptance rate
     */
    double generate_ensembles(const configuration_type &phiStart, size_t amount_of_samples,
                              size_t amount_of_thermalization_steps, bool expand);

    /**
     * @brief Compute the \p observable_function_pointer of the currently loaded \a ensembles vector
     * @param observable_function_pointer Pointer to function, which returns the value of the observable
     *                                    for the given configuration
     * @return vector of observable
     */
    std::vector<double> compute_observable(double (BaseModel<configuration_type>::*observable_function_pointer)(
            const configuration_type &));

    /**
     * @brief Returns the last element of the \a ensembles vector
     * @param default_phi Default configuration to be returned,
     *        if no configurations are present in the \a ensembles vector
     * @return Last element of the \a ensembles vector
     */
    configuration_type get_last_configuration(const configuration_type &default_phi);

    /**
     * @brief Clear the \a ensembles vector
     */
    void clear_ensembles();

    /**
     * @brief Dumps the HMCGenerator, including the currently loaded \a ensembles vector
     *        and the used \a model to \p root
     * @param root Group which gets to hold the model attributes
     *             and the dataset to which the \a ensembles vector gets dumped to
     * @return Dataset to which the \a ensembles vector was dumped to
     */
    HighFive::DataSet dumpToH5(HighFive::Group &root);

    /**
     * @brief Calculates and dumps the observable \a observable_function_pointer to \p root into the \p name dataset
     * @param observable_function_pointer Pointer to function, which returns the value of the observable
     *                                    for the given configuration
     * @param name Name under which the dataset will be stored
     * @param root Group under which to create the dataset \p name
     * @return Dataset to which the observable was dumped to
     */
    HighFive::DataSet dump_observable(double (
    BaseModel<configuration_type>::*observable_function_pointer)(const configuration_type &),
                                      const std::string &name, HighFive::Group &root);

    /**
     * @brief Returns the dataset inside \p root, which holds an ensembles vector
     * @param root Group under which th ensembles vector dataset should reside
     * @return Ensembles vector dataset
     */
    HighFive::DataSet get_dataset(HighFive::Group &root);

private:
    /**
     * @brief Do one HMC step
     * @param phi0 Starting field
     * @return New field (after accept/reject)
     */
    configuration_type do_HMC_step([[maybe_unused]] [[maybe_unused]] const configuration_type &phi0);

    /**
     * @brief Model used in the HMC evolution
     */
    BaseModel<configuration_type> &model;

    /**
     * @brief Amount of molecular dynamic steps used in the integration process
     */
    size_t amount_of_steps;
    /**
     * @brief String to be used as key for the \a amount_of_steps in H5 files
     */
    [[maybe_unused]] static const char *amount_of_steps_name;

    /**
     * @brief Step size to use in the molecular dynamics integration
     */
    double step_size;
    /**
     * @brief String to be used as key for the \a step_size in H5 files
     */
    [[maybe_unused]] static const char *step_size_name;

    /**
     * @brief Random number generator used for the conjugate momenta sampling
     */
    std::default_random_engine &generator;

    /**
     * @brief Integrator used for the molecular dynamics integration
     */
    LeapFrogIntegrator<configuration_type> integrator;

    /**
     * @brief Array of configurations
     */
    std::vector<configuration_type> ensembles;
    /**
     * @brief String to be used as key for the \a ensembles vector in H5 files
     */
    [[maybe_unused]] static const char *ensembles_name;

    /**
     * @brief Amount of accepted configurations
     */
    int accepted_configurations{0};
};

template<class configuration_type> const char *HMCGenerator<configuration_type>::
        amount_of_steps_name{"amount_of_steps"};
template<class configuration_type> const char *HMCGenerator<configuration_type>::step_size_name{"step_size"};
template<class configuration_type> const char *HMCGenerator<configuration_type>::ensembles_name{"ensembles"};


template<class configuration_type>
configuration_type HMCGenerator<configuration_type>::do_HMC_step([[maybe_unused]] const configuration_type &phi0) {
    configuration_type pi(phi0.rows());
    configuration_type phi(phi0);
    std::normal_distribution<double> gauss(0, 1);
    for ([[maybe_unused]] auto &elem: pi) {
        elem = gauss(generator);
    }
    [[maybe_unused]] double H_start = pi.dot(pi) * 0.5 + model.get_action(phi);
    integrator.integrate(amount_of_steps, step_size, phi, pi); //updates in place
    [[maybe_unused]] double H_end = pi.dot(pi) * 0.5 + model.get_action(phi);

    std::uniform_real_distribution<double> uniformRealDistribution(0., 1.);

    if (exp(H_start - H_end) > uniformRealDistribution(generator)) {
        // Accept
        accepted_configurations++;
        return phi;
    } else {
        // Reject
        return phi0;
    }
}

template<class configuration_type>
double
HMCGenerator<configuration_type>::generate_ensembles(const configuration_type &phiStart, size_t amount_of_samples,
                                                     size_t amount_of_thermalization_steps, bool expand) {
    assert(model.check_dimensions(phiStart));
    configuration_type phi(phiStart);

    int startindex{0};
    if (expand) {
        startindex = ensembles.size();
    }

    ensembles.resize(startindex + amount_of_samples);
    for (int i = 0; i < amount_of_thermalization_steps; ++i) {
        phi = do_HMC_step(phi);
    }
    accepted_configurations = 0;
    for (int i = 0; i < amount_of_samples; ++i) {
        phi = do_HMC_step(phi);
        ensembles[startindex + i] = phi;
    }

    double ret{1.};
    if (amount_of_samples) {
        return ret * accepted_configurations / static_cast<double>(amount_of_samples);
    } else {
        return 0.;
    }


}

template<class configuration_type>
HMCGenerator<configuration_type>::HMCGenerator(BaseModel<configuration_type> &model_, size_t amount_of_steps_,
                                               double step_size_, std::default_random_engine &generator_)
        : model{model_}, amount_of_steps{amount_of_steps_}, step_size{step_size_}, integrator{model_},
          generator{generator_} {
}

template<class configuration_type>
HMCGenerator<configuration_type>::HMCGenerator(BaseModel<configuration_type> &model_, HighFive::Group &root,
                                               std::default_random_engine &generator_)
        :model{model_}, amount_of_steps{}, step_size{}, integrator{model_}, generator{generator_} {

    auto ensembles_dataset = root.getDataSet(ensembles_name);
    ensembles_dataset.getAttribute(amount_of_steps_name).read(amount_of_steps);
    ensembles_dataset.getAttribute(step_size_name).read(step_size);

    model.load_ensemble(ensembles, ensembles_dataset);
}

template<class configuration_type>
HMCGenerator<configuration_type>::HMCGenerator(BaseModel<configuration_type> &model_, HighFive::Group &root,
                                               size_t amount_of_steps_,
                                               double step_size_, std::default_random_engine &generator_)
        :model{model_}, amount_of_steps{amount_of_steps_}, step_size{step_size_}, integrator{model_},
         generator{generator_} {
    auto ensembles_dataset = root.getDataSet(ensembles_name);
    model.load_ensemble(ensembles, ensembles_dataset);
}

template<class configuration_type>
configuration_type HMCGenerator<configuration_type>::get_last_configuration(const configuration_type &default_phi) {
    if (ensembles.empty()) {
        return default_phi;
    } else {
        return ensembles.back();
    }
}

template<class configuration_type>
void HMCGenerator<configuration_type>::clear_ensembles() {
    ensembles.clear();
}

template<class configuration_type>
HighFive::DataSet HMCGenerator<configuration_type>::dumpToH5(HighFive::Group &root) {
    model.dumpToH5(root);

    auto ensemble_dataset = model.dump_ensemble(ensembles, root, ensembles_name);
    write_static_size(amount_of_steps,root,amount_of_steps_name);
    write_static_size(step_size,root,step_size_name);
    return ensemble_dataset;
}

template<class configuration_type>
HighFive::DataSet HMCGenerator<configuration_type>::get_dataset(HighFive::Group &root) {
    return root.getDataSet(ensembles_name);
}

template<class configuration_type>
std::vector<double> HMCGenerator<configuration_type>::compute_observable(
        double (BaseModel<configuration_type>::*observable_function_pointer)(const configuration_type &)) {
    std::vector<double> ret(ensembles.size());
    for (int i = 0; i < ensembles.size(); ++i) {
        ret[i] = (model.*observable_function_pointer)(ensembles[i]);
    }
    return ret;
}

template<class configuration_type>
HighFive::DataSet
HMCGenerator<configuration_type>::dump_observable(
        double (BaseModel<configuration_type>::*observable_function_pointer)(const configuration_type &),
        const std::string &name, HighFive::Group &root) {
    std::vector<size_t> offset;
    HighFive::DataSet target_dataset = add_to_expandable_dataset(
            root, name,
            {ensembles.size()},
            offset);

    std::vector<size_t> count{ensembles.size()};
    target_dataset.select(offset, count).write(compute_observable(observable_function_pointer));
    return target_dataset;
}


#endif //BACHELOR_THESIS_HMCGENERATOR_H
