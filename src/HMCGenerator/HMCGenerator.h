/**
 * @file       HMCGenerator.h
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       24.03.22
 * @todo       Add constructor which loads ensembles and/or the model from a file to extend it
 */


#ifndef BACHELOR_THESIS_HMCGENERATOR_H
#define BACHELOR_THESIS_HMCGENERATOR_H

#include <MyTypes.h>
#include <IsingModel.h> //TODO look at vererbung to generalize this
#include <iostream>

class HMCGenerator {
public:
    explicit HMCGenerator(IsingModel &model_)
            : model{model_} {}


private:
    IsingModel &model;
};


#endif //BACHELOR_THESIS_HMCGENERATOR_H
