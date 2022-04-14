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
#include <random>
#include <MyTypes.h>


class Analyzer {
public:
    explicit Analyzer(HighFive::DataSet &dataset_, std::default_random_engine &generator_);

    Analyzer() = delete;

    std::vector<double> auto_correlation(size_t max_t);

    void block_data(size_t block_size);

    void bootstrap_data(size_t amount_of_sample_sets);

private:
    HighFive::DataSet &dataset;
    HighFive::Group group;
    std::vector<double> data;
    std::vector<double> blocked_data;
    std::vector<double> bootstrapped_data;

    void set_mean();

    double mean;
    static const char *mean_name;

    double bootstrap_mean;
    static const char *bootstrap_mean_name;

    double bootstrap_variance;
    static const char *bootstrap_variance_name;

    std::default_random_engine &generator;
};


#endif //BACHELOR_THESIS_ANALYZER_H
