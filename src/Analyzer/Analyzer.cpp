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

Analyzer::Analyzer(HighFive::Group &group_, const std::string &data_name, int start_index, int end_index,
                   std::default_random_engine &generator_)
        : dataset{group_.getDataSet(data_name)}, group{group_}, int_auto_correlation_time{-1},
          mean{0.}, generator{generator_}, bootstrap_mean{0.}, bootstrap_variance{0.} {
    dataset.read(data);
    assert(!data.empty());
    if (end_index > 0) {
        data.erase(data.begin() + end_index, data.end());
    } else {
        end_index = static_cast<int>(data.size());
    }
    if (start_index > 0) {
        data.erase(data.begin(), data.begin() + start_index);
    } else {
        start_index = 0;
    }

    auto_correlation_name = std::string(auto_correlation_base_name).append("_").append(
            std::to_string(start_index)).append("_").append(std::to_string(end_index));
    int_auto_correlation_time_name = std::string(int_auto_correlation_time_base_name).append("_").append(
            std::to_string(start_index)).append("_").append(std::to_string(end_index));
    int_auto_correlation_time_bias_name = std::string(int_auto_correlation_time_bias_base_name).append("_").append(
            std::to_string(start_index)).append("_").append(std::to_string(end_index));
    int_auto_correlation_time_stat_error_name = std::string(int_auto_correlation_time_stat_error_base_name).append(
            "_").append(std::to_string(start_index)).append("_").append(std::to_string(end_index));
    mean_name = std::string(mean_base_name).append("_").append(
            std::to_string(start_index)).append("_").append(std::to_string(end_index));
    bootstrap_mean_name = std::string(bootstrap_mean_base_name).append("_").append(
            std::to_string(start_index)).append("_").append(std::to_string(end_index));
    bootstrap_variance_name = std::string(bootstrap_variance_base_name).append("_").append(
            std::to_string(start_index)).append("_").append(std::to_string(end_index));

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

    for (double &t: ret) {
        t /= start_value;
    }
    int_auto_correlation_time = 0.5;//because we add gamma[0]=1
    int W = 0;
    int S = 1;
    double g_W = 1;
    double tau_W = 0;
    for (int i = 1; i < ret.size(); ++i) {
        if (g_W < 0) {
            break;
        } else {
            W++;
            int_auto_correlation_time += ret[i];
            tau_W = S / log((2 * int_auto_correlation_time + 1) / (2 * int_auto_correlation_time - 1));
            g_W = exp(-W / tau_W) - tau_W / sqrt(W * static_cast<double>(data.size()));
        }
    }
    VectorX vec_temp(ret.size());
    vec_temp = VectorX::Map(&ret[0], ret.size());
    WriteVectorX(vec_temp, group, auto_correlation_name);
    double int_auto_correlation_time_bias{(2 * W + 1) * int_auto_correlation_time / static_cast<double>(data.size())};
    double int_auto_correlation_time_stat_error{2. * sqrt(W + 0.5 - int_auto_correlation_time) *
                                                int_auto_correlation_time / sqrt(static_cast<double>(data.size()))};

    int_auto_correlation_time /= (1 - (2 * W + 1) / static_cast<double>(data.size()));
    std::cout << "i:\t" << int_auto_correlation_time << "\tW:\t" << W << std::endl;
    // this is really the variance!!!
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
    if (block_size <= 0) {
        if (int_auto_correlation_time > 0) {
            block_size = 2 * static_cast<int>(int_auto_correlation_time + 1);
            if (block_size <= 0) {
                std::cerr << "Block size invalid: " << block_size << std::endl;
                return;
            }
        }
    }
    blocked_data.clear();
    if (size_to_use > 0) {
        blocked_data.resize(size_to_use / block_size);
        bootstrap_mean_name = std::string(bootstrap_mean_name).append("_").append(std::to_string(size_to_use));
        bootstrap_variance_name = std::string(bootstrap_variance_name).append("_").append(std::to_string(size_to_use));
    } else {
        blocked_data.resize(data.size() / block_size);
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
            block_data(static_cast<int>(int_auto_correlation_time + 1) / 2, -1);
        } else {
            block_data(42, -1);
        }
    }
    if (blocked_data.empty()) {
        return;
    }
    if (amount_of_sample_sets < 0) {
        amount_of_sample_sets = static_cast<int>(blocked_data.size());

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
