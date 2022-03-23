/**
 * @file       LeapFrogIntegrator.h
 * @brief      Declarations for Leap Frog
 * @author     nico
 * @version    0.0.1
 * @date       23.03.22
 */


#ifndef BACHELOR_THESIS_LEAPFROGINTEGRATOR_H
#define BACHELOR_THESIS_LEAPFROGINTEGRATOR_H

#include <eigen3/Eigen/Dense>
#include <IsingModel.h> //TODO look at vererbung to generalize this

namespace {
    /**
     * @brief Vector type definition
     */
    typedef Eigen::VectorXd VectorX;
}


/**
 * @brief Implementation of the Leap Frog Integrator
 */
class LeapFrogIntegrator {
public:
    /**
     * @brief Constructor of the Leap Frog Integrator
     * @param N Dimension of the vectors which should be integrated
     * @param model Force function to use in the updates of the momenta
     */
    explicit LeapFrogIntegrator(IsingModel &model_) :
            model{model_} {}

    void integrate(size_t amount_of_steps, double step_size, VectorX &phi, VectorX &pi);


private:
    /**
     * @brief Updates the momentum pi and returns the new momentum pi_new
     * @param phi
     * @param pi
     * @param step_size
     * @return pi_new
     */
    void update_pi(VectorX &phi, VectorX &pi, double step_size);

    /**
     * @brief Updates the field phi and returns the new field phi_new
     * @param phi
     * @param pi
     * @param step_size
     * @return phi_new
     */
    static void update_phi(VectorX &phi, VectorX &pi, double step_size);

    IsingModel &model;
};


#endif //BACHELOR_THESIS_LEAPFROGINTEGRATOR_H
