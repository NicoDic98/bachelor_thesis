/**
 * @file       Analyzer.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       12.04.22
 */

#include "Analyzer.h"

#include <utility>

const char *Analyzer::mean_name{"mean"};

Analyzer::Analyzer(HighFive::DataSet &dataset_)
        : dataset{dataset_}, mean{0.} {
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
    return ret;
}

void Analyzer::set_mean() {
    mean = 0.;
    for (auto elem: data) {
        mean += elem;
    }
    mean /= static_cast<double>(data.size());
    if (dataset.hasAttribute(mean_name)) {
        dataset.getAttribute(mean_name).write(mean);
    } else {
        dataset.createAttribute(mean_name, mean);
    }
}

void Analyzer::block_data(size_t block_size) {
    blocked_data.clear();
    blocked_data.resize(data.size() / block_size);

    for (auto &elem: blocked_data) {
        elem = 0.;
    }

    for (int i = 0; i < blocked_data.size() * block_size; ++i) {
        blocked_data[i / block_size] += data[i];
    }

    for (auto &elem: blocked_data) {
        elem /= static_cast<double>(block_size);
    }
}
