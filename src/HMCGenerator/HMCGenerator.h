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
#include <IsingModel.h> //TODO look at vererbung to generalize this
#include <LeapFrogIntegrator.h>
#include <iostream>
#include <random>

/**
 * @brief Implementation of the standard HMC algorithm
 */
class HMCGenerator {
public:
    /**
     * @brief Constructor of the HMC generator
     * @param model_ Model for which to generate ensembles
     * @param amount_of_steps_ Amount of steps to be used in the integration process
     * @param step_size_ Step size in the integration process
     * @param generator_ Random number generator to be used for the HMC
     */
    HMCGenerator(IsingModel &model_, size_t amount_of_steps_, double step_size_,
                 std::default_random_engine &generator_);

    /**
     * @brief Generate amount_of_samples amount of ensembles, starting from phiStart and doing
     *        amount_of_thermalization_steps thermalization steps in advance
     * @param phiStart Starting field
     * @param amount_of_samples Amount of samples to take
     * @param amount_of_thermalization_steps Amount of thermalization steps
     * @return Acceptance rate
     */
    double generate_ensembles(const VectorX &phiStart,
                              size_t amount_of_samples, size_t amount_of_thermalization_steps = 10);

    /**
     * @brief Compute the magnetization of the currently loaded ensemble
     * @return magnetization
     */
    double compute_magnetization();

    /**
     * @brief Returns beta of the used model
     * @return Inverse Temperature
     */
    double get_beta() const { return model.get_beta(); }


private:
    VectorX do_HMC_step(const VectorX &phi0);

    IsingModel &model;
    size_t amount_of_steps;
    double step_size;
    std::default_random_engine &generator; //TODO Maybe make this static?
    LeapFrogIntegrator integrator;
    std::vector<VectorX> ensembles;
    int accepted_configurations{0};
};


#endif //BACHELOR_THESIS_HMCGENERATOR_H
