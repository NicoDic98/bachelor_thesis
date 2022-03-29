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

#include <utility>

enum class InterpolationType {
    Checkerboard, Black_White
};

template<class configuration_type>
class MultiLevelHMCGenerator {
public:
    MultiLevelHMCGenerator(BaseModel<configuration_type> &model_, std::vector<size_t> nu_pre_,
                           std::vector<size_t> nu_post_, size_t gamma_,
                           InterpolationType inter_type_);

private:
    BaseModel<configuration_type> &model;
    std::vector<size_t> nu_pre;
    std::vector<size_t> nu_post;
    InterpolationType inter_type;
};

template<class configuration_type>
MultiLevelHMCGenerator<configuration_type>::MultiLevelHMCGenerator(BaseModel<configuration_type> &model_,
                                                                   std::vector<size_t> nu_pre_,
                                                                   std::vector<size_t> nu_post_, size_t gamma_,
                                                                   InterpolationType inter_type_)
        : model{model_}, nu_pre{std::move(nu_pre_)}, nu_post{std::move(nu_post_)}, inter_type{inter_type_} {

}


#endif //BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
