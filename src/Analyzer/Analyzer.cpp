/**
 * @file       Analyzer.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       12.04.22
 */

#include "Analyzer.h"


const char *Analyzer::auto_correlation_base_name{"auto_correlation"};
const char *Analyzer::int_auto_correlation_time_base_name{"int_auto_correlation_time"};
const char *Analyzer::int_auto_correlation_time_bias_base_name{"int_auto_correlation_time_bias"};
const char *Analyzer::int_auto_correlation_time_stat_error_base_name{"int_auto_correlation_time_stat_error"};
const char *Analyzer::mean_base_name{"mean"};
const char *Analyzer::bootstrap_mean_base_name{"bootstrap_mean"};
const char *Analyzer::bootstrap_variance_base_name{"bootstrap_variance"};

Analyzer::Analyzer(HighFive::Group &group_, const std::string &data_name,
                   std::default_random_engine &generator_)
        : dataset{group_.getDataSet(data_name)}, group{group_}, int_auto_correlation_time{-1},
          mean{0.}, generator{generator_}, bootstrap_mean{0.}, bootstrap_variance{0.} {
    dataset.read(data);
    assert(!data.empty());

    auto_correlation_name = auto_correlation_base_name;
    int_auto_correlation_time_name = int_auto_correlation_time_base_name;
    int_auto_correlation_time_bias_name = int_auto_correlation_time_bias_base_name;
    int_auto_correlation_time_stat_error_name = int_auto_correlation_time_stat_error_base_name;
    mean_name = mean_base_name;
    bootstrap_mean_name = bootstrap_mean_base_name;
    bootstrap_variance_name = bootstrap_variance_base_name;

    set_mean();
}

std::vector<double> Analyzer::auto_correlation(size_t max_t) {
    if (max_t >= data.size()) {
        max_t = data.size() - 1;
    }
    std::vector<double> ret(max_t + 1, 0.);
    for (int t = 0; t < ret.size(); ++t) {
        for (int i = 0; i < data.size() - t; ++i) {
            ret[t] += (data[i] - mean) * (data[i + t] - mean);
        }
        ret[t] /= static_cast<double>(data.size() - t);
    }
    auto start_value = ret[0];
    int_auto_correlation_time = -0.5;//because we add gamma[0]=1
    for (double &t: ret) {
        t /= start_value;
        int_auto_correlation_time += t;
    }
    VectorX vec_temp(ret.size());
    vec_temp = VectorX::Map(&ret[0], ret.size());
    WriteVectorX(vec_temp, group, auto_correlation_name);

    double int_auto_correlation_time_bias{
            -int_auto_correlation_time * exp(-static_cast<double>(max_t) / int_auto_correlation_time)};
    double int_auto_correlation_time_stat_error{4. * (static_cast<double>(max_t) + 0.5 - int_auto_correlation_time) *
                                                int_auto_correlation_time / static_cast<double>(data.size())};
    write_static_size(int_auto_correlation_time, group, int_auto_correlation_time_name);
    write_static_size(int_auto_correlation_time_bias, group, int_auto_correlation_time_bias_name);
    write_static_size(int_auto_correlation_time_stat_error, group, int_auto_correlation_time_stat_error_name);

    return ret;
}

void Analyzer::set_mean() {
    mean = 0.;
    for (auto elem: data) {
        mean += elem;
    }
    mean /= static_cast<double>(data.size());
    write_static_size(mean, group, mean_name);
}

void Analyzer::block_data(int block_size, int size_to_use) {
    if (block_size < 0) {
        if (int_auto_correlation_time > 0) {
            block_size = 2 * static_cast<int>(int_auto_correlation_time);
            if (block_size < 0) {
                std::cerr << "Block size invalid: " << block_size << std::endl;
            }
        }
    }
    blocked_data.clear();
    if (size_to_use > 0) {
        blocked_data.resize(size_to_use / block_size);
        bootstrap_mean_name = std::string(bootstrap_mean_base_name).append(std::to_string(size_to_use));
        bootstrap_variance_name = std::string(bootstrap_variance_base_name).append(std::to_string(size_to_use));
    } else {
        blocked_data.resize(data.size() / block_size);
        bootstrap_mean_name = bootstrap_mean_base_name;
        bootstrap_variance_name = bootstrap_variance_base_name;
    }

    std::fill(blocked_data.begin(), blocked_data.end(), 0.);

    for (int i = 0; i < blocked_data.size() * block_size; ++i) {
        blocked_data[i / block_size] += data[i];
    }

    for (auto &elem: blocked_data) {
        elem /= static_cast<double>(block_size);
    }
}

void Analyzer::bootstrap_data(int amount_of_sample_sets) {
    if (blocked_data.empty()) {
        if (int_auto_correlation_time > 0) {
            block_data(static_cast<int>(int_auto_correlation_time) / 2, -1);
        } else {
            block_data(42, -1);
        }
    }
    if (blocked_data.empty()) {
        exit(-42);
    }
    if (amount_of_sample_sets < 0) {
        amount_of_sample_sets = blocked_data.size();

    }
    bootstrapped_data.clear();
    bootstrapped_data.resize(amount_of_sample_sets);
    std::fill(bootstrapped_data.begin(), bootstrapped_data.end(), 0.);

    std::uniform_int_distribution<size_t> distribution(0, blocked_data.size() - 1);

    bootstrap_mean = 0.;
    for (auto &bootstrap_est: bootstrapped_data) {
        for (int i = 0; i < blocked_data.size(); ++i) {
            bootstrap_est += blocked_data[distribution(generator)];
        }
        bootstrap_est /= static_cast<double>(blocked_data.size());
        bootstrap_mean += bootstrap_est;
    }
    bootstrap_mean /= static_cast<double>(amount_of_sample_sets);

    bootstrap_variance = 0.;
    for (auto bootstrap_est: bootstrapped_data) {
        bootstrap_variance += (bootstrap_est - bootstrap_mean) * (bootstrap_est - bootstrap_mean);
    }
    bootstrap_variance /= static_cast<double>(amount_of_sample_sets - 1);
    write_static_size(bootstrap_mean, group, bootstrap_mean_name);
    write_static_size(bootstrap_variance, group, bootstrap_variance_name);
}
