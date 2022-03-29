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

enum class InterpolationType {
    Checkerboard, Black_White
};

template<class configuration_type>
class MultiLevelHMCGenerator {
public:
    MultiLevelHMCGenerator(BaseModel<configuration_type> &model_, std::vector<size_t> nu_pre_,
                           std::vector<size_t> nu_post_, size_t gamma_,
                           InterpolationType inter_type_,
                           const std::vector<size_t> &amount_of_steps_, const std::vector<double> &step_sizes_,
                           std::default_random_engine &generator_);

private:
    BaseModel<configuration_type> &model;
    std::vector<size_t> nu_pre;
    std::vector<size_t> nu_post;
    InterpolationType inter_type;
    std::default_random_engine &generator;
    std::vector<HMCGenerator<configuration_type>> HMCStack;
    std::vector<std::unique_ptr<BaseModel<configuration_type>>> ModelStack;
};

template<class configuration_type>
MultiLevelHMCGenerator<configuration_type>::MultiLevelHMCGenerator(BaseModel<configuration_type> &model_,
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

    ModelStack.push_back(std::unique_ptr<BaseModel<configuration_type>>(model.get_copy_of_model()));
    HMCStack.push_back(HMCGenerator(model, amount_of_steps_[0], step_sizes_[0], generator));

    for (int i = 1; i < nu_pre.size(); ++i) {
        ModelStack[i - 1]->print_name();
        ModelStack.push_back(
                std::unique_ptr<BaseModel<configuration_type>>(
                        (ModelStack[i - 1])->get_coarser_model(MatrixX(1, 2))));
        HMCStack.push_back(HMCGenerator(*ModelStack[i], amount_of_steps_[i], step_sizes_[i], generator));
    }
    ModelStack[nu_pre.size() - 1]->print_name();

}


#endif //BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
