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

#include <MyTypes.h>
#include <BaseModel.h>
#include <LeapFrogIntegrator.h>
#include <iostream>
#include <random>

/**
 * @brief Implementation of the standard HMC algorithm
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
                 std::default_random_engine &generator_)
            : model{model_}, amount_of_steps{amount_of_steps_}, step_size{step_size_}, integrator{model_},
              generator{generator_} {
    }

    /**
     * @brief Generate amount_of_samples amount of ensembles, starting from phiStart and doing
     *        amount_of_thermalization_steps thermalization steps in advance
     * @param phiStart Starting field
     * @param amount_of_samples Amount of samples to take
     * @param amount_of_thermalization_steps Amount of thermalization steps
     * @return Acceptance rate
     */
    double generate_ensembles(const configuration_type &phiStart,
                              size_t amount_of_samples, size_t amount_of_thermalization_steps = 10) {
        //TODO add expand option
        assert(model.check_dimensions(phiStart));
        configuration_type phi(phiStart);
        ensembles.resize(amount_of_samples);
        for (int i = 0; i < amount_of_thermalization_steps; ++i) {
            phi = do_HMC_step(phi);
        }
        accepted_configurations = 0;
        for (int i = 0; i < amount_of_samples; ++i) {
            phi = do_HMC_step(phi);
            ensembles[i] = phi;
        }

        double ret{1.};
        return ret * accepted_configurations / amount_of_samples;
    }

    /**
     * @brief Compute the magnetization of the currently loaded ensemble
     * @return magnetization
     */
    double compute_magnetization() {
        double ret{0.};
        for (const auto &elem: ensembles) {
            ret += abs(model.get_magnetization(elem));
        }
        return ret / ensembles.size();
    }

    /**
     * @brief Returns beta of the used model
     * @return Inverse Temperature
     */
    double get_beta() const { return model.get_beta(); }


private:
    configuration_type do_HMC_step(const configuration_type &phi0){
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

    BaseModel<configuration_type> &model;
    size_t amount_of_steps;
    double step_size;
    std::default_random_engine &generator; //TODO Maybe make this static?
    LeapFrogIntegrator<configuration_type> integrator;
    std::vector<configuration_type> ensembles;
    int accepted_configurations{0};
};


#endif //BACHELOR_THESIS_HMCGENERATOR_H
