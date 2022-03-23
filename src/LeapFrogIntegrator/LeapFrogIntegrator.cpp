/**
 * @file       LeapFrogIntegrator.cpp
 * @brief      Implementations for Leap Frog
 * @author     nico
 * @version    0.0.1
 * @date       23.03.22
 */


#include "LeapFrogIntegrator.h"

void LeapFrogIntegrator::update_pi(VectorX &phi, VectorX &pi, double step_size) {
    pi += step_size * model.get_force(phi);
}

void LeapFrogIntegrator::update_phi(VectorX &phi, VectorX &pi, double step_size) {
    phi += step_size * pi;
}

void LeapFrogIntegrator::integrate(size_t amount_of_steps, double step_size, VectorX &phi, VectorX &pi) {
    for (size_t _ = 0; _ < amount_of_steps; _++) {
        update_phi(phi,pi,step_size*0.5);
        update_pi(phi,pi,step_size);
        update_phi(phi,pi,step_size*0.5);
    }
}
