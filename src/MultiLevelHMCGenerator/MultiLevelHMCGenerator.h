/**
 * @file       MultiLevelHMCGenerator.h
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       26.03.22
 */


#ifndef BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
#define BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H

#include <BaseModel.h>
#include <HMCGenerator.h>
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
    MultiLevelHMCGenerator(BaseModel<configuration_type> &model_, std::vector<size_t> nu_pre_,
                           std::vector<size_t> nu_post_, size_t gamma_,
                           InterpolationType InterpolationType_,
                           const std::vector<size_t> &amount_of_steps_, const std::vector<double> &step_sizes_,
                           std::default_random_engine &generator_);

    /**
     * @brief Generate \p amount_of_samples amount of ensembles, starting from \p phiStart and doing
     *        \p amount_of_thermalization_steps thermalization steps in advance
     * @param phiStart Starting field
     * @param amount_of_samples Amount of samples to take
     * @param amount_of_thermalization_steps Amount of thermalization steps
     * @return Acceptance rate
     */
    double generate_ensembles(const configuration_type &phiStart,
                              size_t amount_of_samples, size_t amount_of_thermalization_steps = 10);

    configuration_type LevelRecursion(int level, const configuration_type &phi);

private:
    std::vector<size_t> nu_pre;
    std::vector<size_t> nu_post;
    size_t gamma;
    InterpolationType inter_type;
    std::default_random_engine &generator;
    std::vector<HMCGenerator<configuration_type>> HMCStack;
    std::vector<std::unique_ptr<BaseModel<configuration_type>>> ModelStack;
    std::vector<double> AcceptanceRates;
};

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
    assert(nu_pre.size() == nu_post.size());
    assert(nu_pre.size() == amount_of_steps_.size());
    assert(nu_pre.size() == step_sizes_.size());

    AcceptanceRates.resize(nu_pre.size());
    ModelStack.push_back(std::unique_ptr<BaseModel<configuration_type>>(model_.get_copy_of_model()));
    HMCStack.push_back(HMCGenerator(model_, amount_of_steps_[0], step_sizes_[0], generator));

    for (int i = 1; i < nu_pre.size(); ++i) {
        ModelStack.push_back(
                std::unique_ptr<BaseModel<configuration_type>>(
                        (ModelStack[i - 1])->get_coarser_model(inter_type)));
        HMCStack.push_back(HMCGenerator(*ModelStack[i], amount_of_steps_[i], step_sizes_[i], generator));
    }

}

template<class configuration_type>
double MultiLevelHMCGenerator<configuration_type>::generate_ensembles(const configuration_type &phiStart,
                                                                      size_t amount_of_samples,
                                                                      size_t amount_of_thermalization_steps) {
    for (int i = 0; i < amount_of_thermalization_steps; ++i) {
        LevelRecursion(0, phiStart);
    }
    HMCStack[0].clear_ensembles();
    for (auto &elem: AcceptanceRates) {
        elem = 0.;
    }
    for (int i = 0; i < amount_of_samples; ++i) {
        LevelRecursion(0, phiStart);
    }
    return AcceptanceRates[0] / (amount_of_samples * 2);
}

template<class configuration_type>
configuration_type
MultiLevelHMCGenerator<configuration_type>::LevelRecursion(int level, const configuration_type &phi) {
    configuration_type currentField{phi};
    AcceptanceRates[level] += HMCStack[level].generate_ensembles(currentField, nu_pre[level], 0, level == 0);
    currentField = HMCStack[level].get_last_configuration();
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
    return HMCStack[level].get_last_configuration();
}


#endif //BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
