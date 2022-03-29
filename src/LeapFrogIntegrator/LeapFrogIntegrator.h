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
 */
template<class configuration_type>
class LeapFrogIntegrator {
public:
    /**
     * @brief Constructor of the Leap Frog Integrator
     * @param N Dimension of the vectors which should be integrated
     * @param model Force function to use in the updates of the momenta
     */
    explicit LeapFrogIntegrator(BaseModel<configuration_type> &model_) : model{model_} {}

    /**
     * @brief Integrates the field phi with momentum pi along the Hamiltonian equations of motion.
     * @note  Integration happens inplace.
     * @param amount_of_steps Amount of integration steps
     * @param step_size Step size for each integration step
     * @param phi Field
     * @param pi Momentum
     */
    void integrate(size_t amount_of_steps, double step_size, configuration_type &phi, configuration_type &pi){
        for (size_t _ = 0; _ < amount_of_steps; _++) {
            model.update_phi(phi,pi,step_size*0.5);
            model.update_pi(phi,pi,step_size);
            model.update_phi(phi,pi,step_size*0.5);
        }
    }



private:
    BaseModel<configuration_type> &model;
};


#endif //BACHELOR_THESIS_LEAPFROGINTEGRATOR_H
