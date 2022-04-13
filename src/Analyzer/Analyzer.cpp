/**
 * @file       Analyzer.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       12.04.22
 */

#include "Analyzer.h"


const char *Analyzer::mean_name{"mean"};
const char *Analyzer::bootstrap_mean_name{"bootstrap_mean"};
const char *Analyzer::bootstrap_variance_name{"bootstrap_variance"};

Analyzer::Analyzer(HighFive::DataSet &dataset_, std::default_random_engine &generator_)
        : dataset{dataset_}, mean{0.}, generator{generator_},
          bootstrap_mean{0.}, bootstrap_variance{0.} {
    dataset.read(data);
    assert(!data.empty());
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
        ret[t] /= ret[0];
    }
    VectorX vec_temp(ret.size());
    vec_temp = VectorX::Map(&ret[0], ret.size());

    return ret;
}

void Analyzer::set_mean() {
    mean = 0.;
    for (auto elem: data) {
        mean += elem;
    }
    mean /= static_cast<double>(data.size());
    write_static_size(mean, dataset, mean_name);
}

void Analyzer::block_data(size_t block_size) {
    blocked_data.clear();
    blocked_data.resize(data.size() / block_size);
    std::fill(blocked_data.begin(), blocked_data.end(), 0.);

    for (int i = 0; i < blocked_data.size() * block_size; ++i) {
        blocked_data[i / block_size] += data[i];
    }

    for (auto &elem: blocked_data) {
        elem /= static_cast<double>(block_size);
    }
}

void Analyzer::bootstrap_data(size_t amount_of_sample_sets) {
    if (blocked_data.empty()) {
        block_data(42);
    }
    if (blocked_data.empty()) {
        exit(-42);
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
    write_static_size(bootstrap_mean, dataset, bootstrap_mean_name);
    write_static_size(bootstrap_variance, dataset, bootstrap_variance_name);
}
