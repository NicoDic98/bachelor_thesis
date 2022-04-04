/**
 * @file       HMCGenerator.h
 * @brief      Declarations of standard HMC
 * @author     nico
 * @version    0.0.1
 * @date       24.03.22
 * @todo       Add constructor which loads ensembles and/or the model from a file to extend it
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
     * @brief Constructor of the HMC generator
     * @param model_ Model for which to generate ensembles
     * @param amount_of_steps_ Amount of steps to be used in the integration process
     * @param step_size_ Step size in the integration process
     * @param generator_ Random number generator to be used for the HMC
     */
    HMCGenerator(BaseModel<configuration_type> &model_, size_t amount_of_steps_, double step_size_,
                 std::default_random_engine &generator_);

    HMCGenerator(BaseModel<configuration_type> &model_, HighFive::File &file,
                 const std::string &path, std::default_random_engine &generator_);

    HMCGenerator(BaseModel<configuration_type> &model_, HighFive::File &file,
                 const std::string &path, size_t amount_of_steps_, double step_size_,
                 std::default_random_engine &generator_);

    /**
     * @brief Generate amount_of_samples amount of ensembles, starting from phiStart and doing
     *        amount_of_thermalization_steps thermalization steps in advance
     * @param phiStart Starting field
     * @param amount_of_samples Amount of samples to take
     * @param amount_of_thermalization_steps Amount of thermalization steps
     * @return Acceptance rate
     */
    double generate_ensembles(const configuration_type &phiStart, size_t amount_of_samples,
                              size_t amount_of_thermalization_steps, bool expand);

    /**
     * @brief Compute the \p observable_function_pointer of the currently loaded ensemble
     * @return vector of observable
     */
    std::vector<double> compute_observable(double (BaseModel<configuration_type>::*observable_function_pointer)(
            const configuration_type &));

    /**
     * @brief Returns beta of the used model
     * @return Inverse Temperature
     */
    [[maybe_unused]] [[nodiscard]] double get_beta() const { return model.get_beta(); }


    /**
     * @brief Returns the last element of \a ensembles
     * @return Last element of \a ensembles
     */
    configuration_type get_last_configuration(const configuration_type &default_phi);

    /**
     * @brief Clear \a ensembles
     */
    void clear_ensembles();

    void dumpToH5(HighFive::Group &root, std::string sub_name = std::string("ensembles"));

    /**
     * @brief Dump the \p observable_function_pointer of the currently loaded ensemble to \p file at \p path
     */
    void dump_observable(double (
    BaseModel<configuration_type>::*observable_function_pointer)(const configuration_type &),
                         const std::string &path, HighFive::File &file);


private:
    /**
     * @brief Do one HMC step
     * @param phi0 Starting field
     * @return New field (after accept/reject)
     */
    configuration_type do_HMC_step([[maybe_unused]] const configuration_type &phi0);

    /**
     * @brief Model used in the HMC evolution
     */
    BaseModel<configuration_type> &model;

    /**
     * @brief Amount of molecular dynamic steps used in the integration process
     */
    size_t amount_of_steps;
    [[maybe_unused]] static const char *amount_of_steps_name;

    /**
     * @brief Step size to use in the molecular dynamics integration
     */
    double step_size;
    [[maybe_unused]] static const char *step_size_name;

    /**
     * @brief Random number generator used for the conjugate momenta sampling
     */
    std::default_random_engine &generator; //TODO Maybe make this static?

    /**
     * @brief Integrator used for the molecular dynamics integration
     */
    LeapFrogIntegrator<configuration_type> integrator;

    /**
     * @brief Array of configurations
     */
    std::vector<configuration_type> ensembles;

    /**
     * @brief Amount of accepted configurations
     */
    int accepted_configurations{0};
};

template<class configuration_type>
configuration_type HMCGenerator<configuration_type>::do_HMC_step(const configuration_type &phi0) {
    configuration_type pi(phi0.rows());
    configuration_type phi(phi0);
    std::normal_distribution<double> gauss(0, 1);
    for (auto &elem: pi) {
        elem = gauss(generator);
    }
    double H_start = pi.dot(pi) * 0.5 + model.get_action(phi);
    integrator.integrate(amount_of_steps, step_size, phi, pi); //updates in place
    double H_end = pi.dot(pi) * 0.5 + model.get_action(phi);

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
        return ret * accepted_configurations / amount_of_samples;
    } else {
        return 0.;
    }


}

template<class configuration_type> const char *HMCGenerator<configuration_type>::amount_of_steps_name{
        "amount_of_steps"};
template<class configuration_type> const char *HMCGenerator<configuration_type>::step_size_name{"step_size"};

template<class configuration_type>
HMCGenerator<configuration_type>::HMCGenerator(BaseModel<configuration_type> &model_, size_t amount_of_steps_,
                                               double step_size_, std::default_random_engine &generator_)
        : model{model_}, amount_of_steps{amount_of_steps_}, step_size{step_size_}, integrator{model_},
          generator{generator_} {
}

template<class configuration_type>
HMCGenerator<configuration_type>::HMCGenerator(BaseModel<configuration_type> &model_, HighFive::File &file,
                                               const std::string &path,
                                               std::default_random_engine &generator_)
        :model{model_}, integrator{model_}, generator{generator_} {

    amount_of_steps = H5Easy::loadAttribute<size_t>(file, path, amount_of_steps_name);
    step_size = H5Easy::loadAttribute<double>(file, path, step_size_name);

    model.load_ensemble(ensembles, file, path);
}

template<class configuration_type>
HMCGenerator<configuration_type>::HMCGenerator(BaseModel<configuration_type> &model_, HighFive::File &file,
                                               const std::string &path, size_t amount_of_steps_,
                                               double step_size_, std::default_random_engine &generator_)
        :model{model_}, amount_of_steps{amount_of_steps_}, step_size{step_size_}, integrator{model_},
         generator{generator_} {

    model.load_ensemble(ensembles, file, path);
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
void HMCGenerator<configuration_type>::dumpToH5(HighFive::Group &root, std::string sub_name) {
    auto path = root.getPath();
    model.dumpToH5(root);

    if (path.back() != '/') {
        sub_name.insert(0, "/");

    }
    path.append(sub_name);
    auto ensemble_dataset = H5Easy::dump(root.getFile(), path, ensembles);
    ensemble_dataset.createAttribute(amount_of_steps_name, amount_of_steps);
    ensemble_dataset.createAttribute(step_size_name, step_size);

}

template<class configuration_type>
std::vector<double> HMCGenerator<configuration_type>::compute_observable(
        double (BaseModel<configuration_type>::*observable_function_pointer)(const configuration_type &)) {
    std::cout << ensembles.size() << std::endl;
    std::vector<double> ret(ensembles.size());
    for (int i = 0; i < ensembles.size(); ++i) {
        ret[i] = (model.*observable_function_pointer)(ensembles[i]);
    }
    return ret;
}

template<class configuration_type>
void
HMCGenerator<configuration_type>::dump_observable(
        double (BaseModel<configuration_type>::*observable_function_pointer)(const configuration_type &),
        const std::string &path, HighFive::File &file) {
    H5Easy::dump(file, path, compute_observable(observable_function_pointer));
}


#endif //BACHELOR_THESIS_HMCGENERATOR_H
