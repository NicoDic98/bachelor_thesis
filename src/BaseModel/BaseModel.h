/**
 * @file       BaseModel.h
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       27.03.22
 */


#ifndef BACHELOR_THESIS_BASEMODEL_H
#define BACHELOR_THESIS_BASEMODEL_H

template<class configuration_type>
class BaseModel {
public:
    explicit BaseModel(double beta_) : beta{beta_} {}

    virtual ~BaseModel() = default;

    virtual double get_action(const configuration_type&) = 0;

    /**
     * @brief Returns beta
     * @return Inverse temperature
     */
    virtual double get_beta() {
        return beta;
    }

    /**
     * @brief Sets the value of beta to new_beta
     * @param new_beta Inverse temperature
     */
    virtual void set_beta(double new_beta) {
        beta = new_beta;
    }

    virtual double get_magnetization(const configuration_type&) = 0;


    virtual bool check_dimensions(const configuration_type&) = 0;

    /**
     * @brief Updates the momentum pi and returns the new momentum pi_new
     * @param phi
     * @param pi
     * @param step_size
     * @return pi_new
     */
    virtual void update_pi(configuration_type &phi, configuration_type &pi, double step_size){
        pi += step_size * get_force(phi);
    }

    /**
     * @brief Updates the field phi and returns the new field phi_new
     * @param phi
     * @param pi
     * @param step_size
     * @return phi_new
     */
    virtual void update_phi(configuration_type &phi, configuration_type &pi, double step_size){
        phi += step_size * pi;
    }

    virtual configuration_type get_dof_sample(){
        return configuration_type();
    }

protected:

    double beta;

    virtual configuration_type get_force(const configuration_type&) = 0;

};


#endif //BACHELOR_THESIS_BASEMODEL_H
