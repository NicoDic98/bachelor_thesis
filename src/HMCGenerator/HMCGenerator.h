/**
 * @file       HMCGenerator.h
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       24.03.22
 */


#ifndef BACHELOR_THESIS_HMCGENERATOR_H
#define BACHELOR_THESIS_HMCGENERATOR_H

#include <IsingModel.h> //TODO look at vererbung to generalize this
#include <iostream>

class HMCGenerator {
public:
    HMCGenerator(IsingModel &model_)
            : model{model_} {}

private:
    IsingModel &model;
};


#endif //BACHELOR_THESIS_HMCGENERATOR_H
