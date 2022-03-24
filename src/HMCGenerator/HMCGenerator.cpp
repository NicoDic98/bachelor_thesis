/**
 * @file       HMCGenerator.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       24.03.22
 */


#include "HMCGenerator.h"


std::default_random_engine HMCGenerator::generator{42};

HMCGenerator::HMCGenerator(IsingModel &model_, size_t amount_of_steps_, double step_size_)
        : model{model_}, amount_of_steps{amount_of_steps_}, step_size{step_size_}, integrator{model_} {
}

VectorX HMCGenerator::do_HMC_step(const VectorX &phi0) {
    VectorX pi(phi0.rows());
    VectorX phi(phi0);
    std::normal_distribution<double> gauss(0, 1);
    for (auto &elem: pi) {
        elem = gauss(generator);
    }
    double H_start = pi.dot(pi) * 0.5 + model.get_action(phi);
    integrator.integrate(amount_of_steps, step_size, phi, pi); //updates in place
    double H_end = pi.dot(pi) * 0.5 + model.get_action(phi);

    std::uniform_real_distribution<double> uniformRealDistribution(0., 1.);

    if ((H_start - H_end) > uniformRealDistribution(generator)) {
        // Accept
        acceptance.push_back(true);
        return phi;
    } else {
        acceptance.push_back(false);
        // Reject
        return phi0;
    }
}

double HMCGenerator::generate_ensembles(const VectorX &phiStart,
                                        size_t amount_of_samples, size_t amount_of_thermalization_steps) {
    //TODO add expand option
    VectorX phi(phiStart);
    for (int i = 0; i < amount_of_thermalization_steps; ++i) {
        phi = do_HMC_step(phi);
    }

    for (int i = 0; i < amount_of_thermalization_steps; ++i) {
        phi = do_HMC_step(phi);
        ensembles.push_back(phi);
    }

    double ret{0.};
    for (auto elem: acceptance) {
        ret += elem;
    }
    ret = ret / acceptance.size();
    return ret;
}
