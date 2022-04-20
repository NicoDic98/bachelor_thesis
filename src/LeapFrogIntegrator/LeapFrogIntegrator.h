/**
 * @file       LeapFrogIntegrator.h
 * @brief      Declarations for Leap Frog
 * @author     nico
 * @version    0.0.1
 * @date       23.03.22
 */


#ifndef BACHELOR_THESIS_LEAPFROGINTEGRATOR_H
#define BACHELOR_THESIS_LEAPFROGINTEGRATOR_H

#include <BaseModel.h>

/**
 * @brief Implementation of the Leap Frog Integrator
 * @tparam configuration_type Datatype, which is used for the configurations of the model
 */
template<class configuration_type>
class LeapFrogIntegrator {
public:
    /**
     * @brief Constructor of the Leap Frog Integrator
     * @param model_ Model to be used for integration
     */
    explicit LeapFrogIntegrator(BaseModel<configuration_type> &model_);

    /**
     * @brief Integrates the field \p phi with momentum \p pi along the Hamiltonian equations of motion.
     * @note  Integration happens inplace.
     * @param amount_of_steps Amount of integration steps
     * @param step_size Step size for each integration step
     * @param phi Field
     * @param pi Momentum
     */
    void integrate(size_t amount_of_steps, double step_size, configuration_type &phi, configuration_type &pi);


private:
    /**
     * @brief Model used in the integration
     */
    BaseModel<configuration_type> &model;
};

template<class configuration_type>
LeapFrogIntegrator<configuration_type>::LeapFrogIntegrator(BaseModel<configuration_type> &model_) : model{model_} {}

template<class configuration_type>
void LeapFrogIntegrator<configuration_type>::integrate(size_t amount_of_steps, double step_size,
                                                       configuration_type &phi, configuration_type &pi) {
    for (size_t _ = 0; _ < amount_of_steps; _++) {
        model.update_phi(phi, pi, step_size * 0.5);
        model.update_pi(phi, pi, step_size);
        model.update_phi(phi, pi, step_size * 0.5);
    }
}


#endif //BACHELOR_THESIS_LEAPFROGINTEGRATOR_H
