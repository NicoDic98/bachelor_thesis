/**
 * @file       HMCGenerator.h
 * @brief
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

class HMCGenerator {
public:
    HMCGenerator(IsingModel &model_, size_t amount_of_steps_, double step_size_,
                 std::default_random_engine &generator_);

    double generate_ensembles(const VectorX &phiStart,
                              size_t amount_of_samples, size_t amount_of_thermalization_steps = 10);

    double compute_magnetization();

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
