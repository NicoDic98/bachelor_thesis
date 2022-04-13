/**
 * @file       MultiLevelHMCGenerator.h
 * @brief      Declarations of Multi level HMC
 * @author     nico
 * @version    0.0.2
 * @date       26.03.22
 */


#ifndef BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
#define BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H

#include <BaseModel.h>
#include <HMCGenerator.h>
#include <Analyzer.h>
#include <random>
#include <iostream>

#include <utility>
#include <memory>

/**
 * @brief Template of the Multilevel HMC algorithm
 * @tparam configuration_type Datatype, which is used for the configurations of the model
 */
template<class configuration_type>
class MultiLevelHMCGenerator {
public:
    /**
     * @brief Constructor of the Multi Level HMC generator
     * @param model_ Model for which to generate ensembles
     * @param nu_pre_ Amount of pre coarsening steps to take at each level
     * @param nu_post_ Amount of post coarsening steps to take at each level
     * @param gamma_ Amount of repetitions at each level (determines if a 'V' or 'W' or ... cycle is performed)
     * @param InterpolationType_ Interpolation type used to generate the coarser levels
     * @param amount_of_steps_ Amount of steps to be used in the integration process for each level
     * @param step_sizes_ Step size in the integration process for each level
     * @param generator_ Random number generator to be used for the HMC process
     */
    MultiLevelHMCGenerator(BaseModel<configuration_type> &model_, std::vector<size_t> nu_pre_,
                           std::vector<size_t> nu_post_, size_t gamma_,
                           InterpolationType InterpolationType_,
                           const std::vector<size_t> &amount_of_steps_, const std::vector<double> &step_sizes_,
                           std::default_random_engine &generator_);

    /**
     * @brief Loads the MultiLevelHMCGenerator from \ p file
     * @param model_ Model for which to generate ensembles
     * @param file File which has the levels as groups, which hold the level attributes
     * @param generator_ Random number generator to be used for the HMC process
     */
    MultiLevelHMCGenerator(BaseModel<configuration_type> &model_,
                           HighFive::File &file,
                           std::default_random_engine &generator_);

    /**
     * @brief Generate \p amount_of_samples amount of ensembles, starting from \p phiStart and doing
     *        \p amount_of_thermalization_steps thermalization steps in advance
     * @param phiStart Starting field
     * @param amount_of_samples Amount of samples to take
     * @param amount_of_thermalization_steps Amount of thermalization steps
     * @return Acceptance rates
     */
    std::vector<double> generate_ensembles(const configuration_type &phiStart,
                                           size_t amount_of_samples, size_t amount_of_thermalization_steps = 10);

    /**
     * @brief Calculates and dumps the observable \a observable_function_pointer for each level to \p file
     *        into subgroups for each level intern into the \p name dataset
     * @param observable_function_pointer Pointer to function, which returns the value of the observable
     *                                    for the given configuration
     * @param name Name under which the dataset will be stored
     * @param file File to dump to
     */
    void dump_observable(double (
    BaseModel<configuration_type>::*observable_function_pointer)(const configuration_type &),
                         const std::string &name, HighFive::File &file);

    /**
     * @brief Dumps the MultiLevelHMCGenerator, including all sub models and sub HMCGenerators
     *        into different subgroups with corresponding level parameters
     * @param file File to dump to
     */
    void dumpToH5(HighFive::File &file);

    /**
     * @brief Propagates the change of attributes from the fine level model to all coarser levels
     */
    [[maybe_unused]] void propagate_update();


private:
    /**
     * @brief Recursion to go to coarser levels
     * @param level Current level id (finest=0, coarsest=\c nu_pre.size()-1 )
     * @param phi Starting configuration/field
     * @return Updated configuration/field
     */
    configuration_type LevelRecursion(int level, const configuration_type &phi);

    /**
     * @brief Prefix for the level groups in H5 files
     */
    [[maybe_unused]] static const char *level_name;

    /**
     * @brief Name for the measurements groups in H5 files
     */
    [[maybe_unused]] static const char *measurements_name;

    /**
     * @brief Amount of pre coarsening steps to take at each level
     */
    std::vector<size_t> nu_pre;
    /**
     * @brief String to be used as key for the \a nu_pre in H5 files
     */
    [[maybe_unused]] static const char *nu_pre_name;

    /**
     * @brief Amount of post coarsening steps to take at each level
     */
    std::vector<size_t> nu_post;
    /**
     * @brief String to be used as key for the \a nu_post in H5 files
     */
    [[maybe_unused]] static const char *nu_post_name;

    /**
     * @brief Amount of repetitions at each level (determines if a 'V' or 'W' or ... cycle is performed)
     */
    size_t gamma;
    /**
     * @brief String to be used as key for the \a gamma in H5 files
     */
    [[maybe_unused]] static const char *gamma_name;

    /**
     * @brief Interpolation type used to generate the coarser levels
     */
    InterpolationType inter_type;
    /**
     * @brief String to be used as key for the \a inter_type in H5 files
     */
    [[maybe_unused]] static const char *inter_type_name;

    /**
     * @brief Random number generator to be used for the HMC process
     */
    std::default_random_engine &generator;

    /**
     * @brief Stores the HMC generator for each level
     */
    std::vector<HMCGenerator<configuration_type>> HMCStack;

    /**
     * @brief Stores the model for each level
     */
    std::vector<std::unique_ptr<BaseModel<configuration_type>>> ModelStack;

    /**
     * @brief Saves the Acceptance rate at each level
     */
    std::vector<double> AcceptanceRates;
    /**
     * @brief String to be used as key for the \a AcceptanceRate in H5 files
     */
    [[maybe_unused]] static const char *AcceptanceRate_name;
};

template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::level_name{"level"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::
        measurements_name{"measurements"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::nu_pre_name{"nu_pre"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::nu_post_name{"nu_post"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::gamma_name{"gamma"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::
        inter_type_name{"inter_type"};
template<class configuration_type> const char *MultiLevelHMCGenerator<configuration_type>::
        AcceptanceRate_name{"AcceptanceRate"};


template<class configuration_type>
MultiLevelHMCGenerator<configuration_type>::MultiLevelHMCGenerator(BaseModel<configuration_type> &model_,
                                                                   std::vector<size_t> nu_pre_,
                                                                   std::vector<size_t> nu_post_,
                                                                   size_t gamma_,
                                                                   InterpolationType InterpolationType_,
                                                                   const std::vector<size_t> &amount_of_steps_,
                                                                   const std::vector<double> &step_sizes_,
                                                                   std::default_random_engine &generator_)
        : nu_pre{std::move(nu_pre_)}, nu_post{std::move(nu_post_)}, gamma{gamma_}, inter_type{InterpolationType_},
          generator{generator_}, AcceptanceRates{} {
    //TODO: add auto sizing
    assert(gamma > 0);
    assert(nu_pre.size() == nu_post.size());
    assert(nu_pre.size() == amount_of_steps_.size());
    assert(nu_pre.size() == step_sizes_.size());


    AcceptanceRates.resize(nu_pre.size());
    ModelStack.push_back(std::unique_ptr<BaseModel<configuration_type>>(model_.get_copy_of_model()));
    HMCStack.push_back(HMCGenerator(*ModelStack[0], amount_of_steps_[0], step_sizes_[0], generator));

    for (int i = 1; i < nu_pre.size(); ++i) {
        assert(nu_pre[i] + nu_post[i] > 0);
        ModelStack.push_back(
                std::unique_ptr<BaseModel<configuration_type>>(
                        (ModelStack[i - 1])->get_coarser_model(inter_type)));
        HMCStack.push_back(HMCGenerator(*ModelStack[i], amount_of_steps_[i], step_sizes_[i], generator));
    }

}

template<class configuration_type>
MultiLevelHMCGenerator<configuration_type>::MultiLevelHMCGenerator(BaseModel<configuration_type> &model_,
                                                                   HighFive::File &file,
                                                                   std::default_random_engine &generator_)
        : gamma{}, generator{generator_}, AcceptanceRates{} {

    std::string current_level_name{level_name};
    current_level_name.append(std::to_string(0));
    assert(file.exist(current_level_name));

    HighFive::Group current_level = file.getGroup(current_level_name);

    ModelStack.push_back(std::unique_ptr<BaseModel<configuration_type>>(model_.get_copy_of_model()));
    HMCStack.push_back(HMCGenerator(*ModelStack[0], current_level, generator));
    HighFive::DataSet current_level_dataset = HMCStack[0].get_dataset(current_level);

    size_t buffer;
    current_level_dataset.getAttribute(nu_pre_name).read(buffer);
    nu_pre.push_back(buffer);
    current_level_dataset.getAttribute(nu_post_name).read(buffer);
    nu_post.push_back(buffer);
    current_level_dataset.getAttribute(gamma_name).read(gamma);
    double dBuffer;
    current_level_dataset.getAttribute(AcceptanceRate_name).read(dBuffer);
    AcceptanceRates.push_back(dBuffer);
    assert(nu_pre[0] + nu_post[0] > 0);

    for (int i = 1; i < 1000; ++i) {
        current_level_name = level_name;
        current_level_name.append(std::to_string(i));
        if (!file.exist(current_level_name)) {
            break;
        }

        current_level = file.getGroup(current_level_name);


        ModelStack.push_back(
                std::unique_ptr<BaseModel<configuration_type>>(
                        (ModelStack[i - 1])->get_model_at(current_level)));
        HMCStack.push_back(HMCGenerator(*ModelStack[i], current_level, generator));
        current_level_dataset = HMCStack[i].get_dataset(current_level);

        current_level_dataset.getAttribute(nu_pre_name).read(buffer);
        nu_pre.push_back(buffer);
        current_level_dataset.getAttribute(nu_post_name).read(buffer);
        nu_post.push_back(buffer);
        current_level_dataset.getAttribute(gamma_name).read(gamma);
        int iBuffer;
        current_level_dataset.getAttribute(inter_type_name).read(iBuffer);
        inter_type = static_cast<InterpolationType>(iBuffer);
        current_level_dataset.getAttribute(AcceptanceRate_name).read(dBuffer);
        AcceptanceRates.push_back(dBuffer);
        assert(nu_pre[i] + nu_post[i] > 0);
    }
}

template<class configuration_type>
std::vector<double> MultiLevelHMCGenerator<configuration_type>::generate_ensembles(const configuration_type &phiStart,
                                                                                   size_t amount_of_samples,
                                                                                   size_t amount_of_thermalization_steps) {
    configuration_type phi(phiStart);
    for (int i = 0; i < amount_of_thermalization_steps; ++i) {
        phi = LevelRecursion(0, phi);
    }
    HMCStack[0].clear_ensembles();
    for (auto &elem: AcceptanceRates) {
        elem = 0.;
    }
    for (int i = 0; i < amount_of_samples; ++i) {
        phi = LevelRecursion(0, phi);
    }

    for (int i = 0; i < AcceptanceRates.size(); ++i) {
        AcceptanceRates[i] = AcceptanceRates[i] / (amount_of_samples * (nu_pre[i] + nu_post[i]) * int_pow(gamma, i));
    }
    return AcceptanceRates;
}

template<class configuration_type>
configuration_type
MultiLevelHMCGenerator<configuration_type>::LevelRecursion(int level, const configuration_type &phi) {
    configuration_type currentField{phi};
    AcceptanceRates[level] += HMCStack[level].generate_ensembles(currentField, nu_pre[level], 0, level == 0);
    currentField = HMCStack[level].get_last_configuration(currentField);
    //std::cout << "Level:\t" << level << std::endl;

    if (level < (nu_pre.size() - 1)) {
        ModelStack[level + 1]->update_fields(currentField);
        configuration_type CoarseCorrections = ModelStack[level + 1]->get_empty_field();
        for (int i = 0; i < gamma; ++i) {
            CoarseCorrections = LevelRecursion(level + 1, CoarseCorrections);
        }
        ModelStack[level + 1]->interpolate(CoarseCorrections, currentField);
    }
    AcceptanceRates[level] += HMCStack[level].generate_ensembles(currentField, nu_post[level], 0, level == 0);
    return HMCStack[level].get_last_configuration(currentField);
}

template<class configuration_type>
[[maybe_unused]] void MultiLevelHMCGenerator<configuration_type>::propagate_update() {
    //Start from 0 because this way also level 0 updates its parameters from the root model
    for (int i = 0; i < ModelStack.size(); ++i) {
        ModelStack[i]->pull_attributes_from_root();
    }
}


template<class configuration_type>
void MultiLevelHMCGenerator<configuration_type>::dumpToH5(HighFive::File &file) {
    for (int i = 0; i < HMCStack.size(); ++i) {
        HighFive::Group current_level = file.getGroup(file.getPath());

        std::string current_level_name{level_name};
        current_level_name.append(std::to_string(i));

        if (file.exist(current_level_name)) {
            current_level = file.getGroup(current_level_name);
        } else {
            current_level = file.createGroup(current_level_name);
        }
        auto ensemble_dataset = HMCStack[i].dumpToH5(current_level);
        write_static_size(nu_pre[i], ensemble_dataset, nu_pre_name);
        write_static_size(nu_post[i], ensemble_dataset, nu_post_name);
        write_static_size(gamma, ensemble_dataset, gamma_name);
        write_static_size(static_cast<int>(inter_type), ensemble_dataset, inter_type_name);
        write_static_size(AcceptanceRates[i], ensemble_dataset, AcceptanceRate_name);
    }
}

template<class configuration_type>
void MultiLevelHMCGenerator<configuration_type>::dump_observable(
        double (BaseModel<configuration_type>::*observable_function_pointer)(const configuration_type &),
        const std::string &name, HighFive::File &file) {
    for (int i = 0; i < HMCStack.size(); ++i) {
        HighFive::Group current_level = file.getGroup(file.getPath());
        HighFive::Group measurements = current_level;

        std::string current_level_name{level_name};
        current_level_name.append(std::to_string(i));

        if (file.exist(current_level_name)) {
            current_level = file.getGroup(current_level_name);
        } else {
            current_level = file.createGroup(current_level_name);
        }
        if (current_level.exist(measurements_name)) {
            measurements = current_level.getGroup(measurements_name);
        } else {
            measurements = current_level.createGroup(measurements_name);
        }
        HighFive::DataSet dataset = HMCStack[i].dump_observable(observable_function_pointer, name, measurements);
        Analyzer a(dataset, generator);
        if(i==0){
            a.block_data(16);
            a.bootstrap_data(200);
        }
    }
}


#endif //BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
