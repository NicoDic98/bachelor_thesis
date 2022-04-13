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
    explicit Analyzer(HighFive::DataSet &dataset_);

    Analyzer() = delete;

    std::vector<double> auto_correlation(size_t max_t);

    void block_data(size_t block_size);

private:
    HighFive::DataSet &dataset;
    std::vector<double> data;
    std::vector<double> blocked_data;

    void set_mean();

    double mean;
    static const char *mean_name;
};


#endif //BACHELOR_THESIS_ANALYZER_H
