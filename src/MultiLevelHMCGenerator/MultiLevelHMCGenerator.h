/**
 * @file       MultiLevelHMCGenerator.h
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       26.03.22
 */


#ifndef BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
#define BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H

#include <HMCGenerator.h>

enum class InterpolationType {
    Checkerboard, Black_White
};

class MultiLevelHMCGenerator {
    MultiLevelHMCGenerator(std::vector<size_t> nu_pre_, std::vector<size_t> nu_post_, size_t gamma_,
                           InterpolationType inter_type_);
};


#endif //BACHELOR_THESIS_MULTILEVELHMCGENERATOR_H
