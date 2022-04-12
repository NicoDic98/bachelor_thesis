/**
 * @file       Analyzer.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       12.04.22
 */


#include "Analyzer.h"

#include <utility>

Analyzer::Analyzer(std::vector<double> data_)
        : data{std::move(data_)}, mean{0.} {
    assert(!data.empty());
    set_mean();
}

Analyzer::Analyzer(HighFive::DataSet &dataset_)
        : mean{0.} {
    dataset_.read(data);
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
}
