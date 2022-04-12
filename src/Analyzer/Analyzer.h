/**
 * @file       Analyzer.h
 * @brief      Provides analyzing functionality
 * @author     nico
 * @version    0.0.1
 * @date       12.04.22
 */


#ifndef BACHELOR_THESIS_ANALYZER_H
#define BACHELOR_THESIS_ANALYZER_H

#include <vector>
#include <highfive/H5File.hpp>


class Analyzer {
public:
    explicit Analyzer(std::vector<double> data_);

    explicit Analyzer(HighFive::DataSet& dataset_);

    std::vector<double> auto_correlation(size_t max_t);

private:
    std::vector<double> data;
    void set_mean();
    double mean;
};


#endif //BACHELOR_THESIS_ANALYZER_H
