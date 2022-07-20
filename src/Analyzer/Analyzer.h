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

/**
 * @brief Base object for analyzes
 */
class Analyzer {
public:
    /**
     * @brief Standard constructor for Analyzer
     * @param group_ Group under which the dataset with name \p data_name lies
     * @param data_name Name of dataset
     * @param start_index Index to be used as starting point for the analysis
     * @param end_index Index to be used as ending point for the analysis
     * @param generator_ Random number generator to be used for the bootstrap
     */
    explicit Analyzer(HighFive::Group &group_, const std::string &data_name, int start_index,int end_index,
                      std::default_random_engine &generator_);

    /**
     * @brief Remove default constructor
     */
    Analyzer() = delete;

    /**
     * @brief Calculates the auto correlation up to \p max_t and writes it to \a group as dataset under the name \a auto_correlation_base_name
     * @param max_t Maximum Marcov chain time, up to which the auto correlation gets calculated
     * @return Vector containing the values of the auto correlation
     */
    std::vector<double> auto_correlation(size_t max_t);

    /**
     * @brief Block the data into blocks of size \p block_size
     * @param block_size Block size
     * @param size_to_use Part of internal array to be used for the blocking
     */
    void block_data(int block_size, int size_to_use);

    /**
     * @brief Calculates the bootstrap mean and variance on the previously blocked data
     * @param amount_of_sample_sets Amount of boostrap samples
     */
    void bootstrap_data(int amount_of_sample_sets);

private:
    /**
     * @brief Dataset
     */
    HighFive::DataSet dataset;

    /**
     * @brief Group to which the analyses get written to
     */
    HighFive::Group group;

    /**
     * @brief Vector containing the loaded data
     */
    std::vector<double> data;

    /**
     * @brief Vector containing the blocked data
     */
    std::vector<double> blocked_data;

    /**
     * @brief Vector containing the bootstrap samples
     */
    std::vector<double> bootstrapped_data;

    /**
     * @brief String to be used as key for the auto_correlation in H5 files
     */
    static const char *auto_correlation_base_name;

    /**
     * @brief String to be used as key for the auto_correlation in H5 files
     */
    std::string auto_correlation_name;

    /**
     * @brief Mean of the dataset
     */
    double int_auto_correlation_time;

    /**
     * @brief String to be used as key for the int_auto_correlation_time in H5 files
     */
    static const char *int_auto_correlation_time_base_name;

    /**
     * @brief String to be used as key for the auto_correlation in H5 files
     */
    std::string int_auto_correlation_time_name;

    /**
     * @brief String to be used as key for the int_auto_correlation_time_bias in H5 files
     */
    static const char *int_auto_correlation_time_bias_base_name;

    /**
     * @brief String to be used as key for the int_auto_correlation_time_bias in H5 files
     */
    std::string int_auto_correlation_time_bias_name;

    /**
     * @brief String to be used as key for the int_auto_correlation_time_stat_error in H5 files
     */
    static const char *int_auto_correlation_time_stat_error_base_name;

    /**
     * @brief String to be used as key for the int_auto_correlation_time_stat_error in H5 files
     */
    std::string int_auto_correlation_time_stat_error_name;

    /**
     * @brief Calculates, sets and saves the mean of the dataset
     */
    void set_mean();

    /**
     * @brief Mean of the dataset
     */
    double mean;
    /**
     * @brief String to be used as key for \a mean in H5 files
     */
    static const char *mean_base_name;
    /**
     * @brief String to be used as key for \a mean in H5 files
     */
    std::string mean_name;

    /**
     * @brief Bootstrap mean of the dataset
     */
    double bootstrap_mean;
    /**
     * @brief String to be used as key for \a bootstrap_mean in H5 files
     */
    static const char *bootstrap_mean_base_name;
    /**
     * @brief String to be used as key for \a bootstrap_mean in H5 files
     */
    std::string bootstrap_mean_name;

    /**
     * @brief Bootstrap variance of the dataset
     */
    double bootstrap_variance;
    /**
     * @brief String to be used as key for \a bootstrap_variance in H5 files
     */
    static const char *bootstrap_variance_base_name;
    /**
     * @brief String to be used as key for \a bootstrap_variance in H5 files
     */
    std::string bootstrap_variance_name;

    /**
     * @brief Random number generator to be used for the bootstrap
     */
    std::default_random_engine &generator;
};


#endif //BACHELOR_THESIS_ANALYZER_H
