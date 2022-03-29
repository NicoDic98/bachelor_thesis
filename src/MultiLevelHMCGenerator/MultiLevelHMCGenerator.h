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

enum class InterpolationType {
    Checkerboard, Black_White
};

template<class configuration_type, class model_type>
class MultiLevelHMCGenerator {
public:
    MultiLevelHMCGenerator(model_type &model_, std::vector<size_t> nu_pre_,
                           std::vector<size_t> nu_post_, size_t gamma_,
                           InterpolationType inter_type_,
                           const std::vector<size_t> &amount_of_steps_, const std::vector<double> &step_sizes_,
                           std::default_random_engine &generator_);

private:
    model_type &model;
    std::vector<size_t> nu_pre;
    std::vector<size_t> nu_post;
    InterpolationType inter_type;
    std::default_random_engine &generator;
    std::vector<HMCGenerator<configuration_type>> HMCStack;
    std::vector<model_type> ModelStack;
};

template<class configuration_type, class model_type>
MultiLevelHMCGenerator<configuration_type, model_type>::MultiLevelHMCGenerator(model_type &model_,
                                                                   std::vector<size_t> nu_pre_,
                                                                   std::vector<size_t> nu_post_,
                                                                   size_t gamma_,
                                                                   InterpolationType inter_type_,
                                                                   const std::vector<size_t> &amount_of_steps_,
                                                                   const std::vector<double> &step_sizes_,
                                                                   std::default_random_engine &generator_)
        : model{model_}, nu_pre{std::move(nu_pre_)}, nu_post{std::move(nu_post_)}, inter_type{inter_type_},
          generator{generator_} {
    //TODO: add auto sizing
    assert(nu_pre.size() == nu_post.size());
    assert(nu_pre.size() == amount_of_steps_.size());
    assert(nu_pre.size() == step_sizes_.size());
    HMCStack.push_back(HMCGenerator(model, amount_of_steps_[0], step_sizes_[0], generator));
    decltype(model) a(model);
    a.print_name();
    ModelStack.push_back(model);
    for (int i = 1; i < nu_pre.size(); ++i) {
        //ModelStack.push_back(decltype(model_)());
        HMCStack.push_back(HMCGenerator(model, amount_of_steps_[i], step_sizes_[i], generator));
    }
}


#endif //BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
