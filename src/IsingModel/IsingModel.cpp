/**
 * @file       IsingModel.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       23.03.22
 */


#include <iostream>
#include <utility>
#include "IsingModel.h"


IsingModel::IsingModel(double beta_, VectorX h_, VectorX eta_, double offset_, int dimension_,
                       int neighbour_extent_, int grid_size_)
        : beta{beta_}, sqrt_beta{sqrt(beta_)}, h{std::move(h_)}, eta{std::move(eta_)}, offset{offset_},
          k_sym(int_pow(grid_size_, dimension_), int_pow(grid_size_, dimension_)),
          k_rec(int_pow(grid_size_, dimension_), int_pow(grid_size_, dimension_)) {
    fill_connectivity_matrix(dimension_, neighbour_extent_, grid_size_);
    add_offset_to_connectivity_matrix();
    k_rec = k_sym;
}

IsingModel::IsingModel(double beta_, VectorX h_, VectorX eta_, double offset_, MatrixX k_sym_, MatrixX k_rec_)
        : beta{beta_}, sqrt_beta{sqrt(beta_)}, h{std::move(h_)}, eta{std::move(eta_)}, offset{offset_},
          k_sym{std::move(k_sym_)},
          k_rec{std::move(k_rec_)} {


}

double IsingModel::get_action(const VectorX &phi) {


    return 0;
}

VectorX IsingModel::get_force(const VectorX &phi) {
    VectorX phi_tilde{k_rec * phi};
    VectorX temp{eta + sqrt_beta * phi_tilde};
    for (auto &elem: temp) {
        elem = tanh(elem);
    }
    return -phi_tilde + sqrt_beta * h + sqrt_beta * k_rec.transpose() * temp;
}

double IsingModel::get_magnetization(const VectorX &phi) {
    return 0;
}

void IsingModel::fill_connectivity_matrix(int dimension, int neighbour_extent, int grid_size) {
    k_sym.setZero();

    std::vector<long> index_offsets;
    long size = k_sym.rows();
    generate_index_offsets(index_offsets, size, 0, 0, 0);

    for (auto elem: index_offsets) {
        std::cout << elem << '\n';
    }

    for (long m = 0; m < size; ++m) {
        for (auto index_offset: index_offsets) {
            k_sym(m, (m + index_offset) % size) += 1;
        }

    }
}

void IsingModel::add_offset_to_connectivity_matrix() {
    MatrixX OffsetMatrix(k_sym.rows(), k_sym.cols());
    OffsetMatrix.setIdentity();
    k_sym += offset * OffsetMatrix;
}

void IsingModel::print_connectivity_matrix() {
    std::cout << k_sym;
}

void
IsingModel::generate_index_offsets(std::vector<long> &index_offsets, long size, int dimension, int neighbour_extent,
                                   int grid_size) {
    long base_offset{1};
    for (int i = 0; i < dimension; ++i) {
        for (long j = 1; j <= neighbour_extent; ++j) {
            index_offsets.push_back(j * base_offset);
            index_offsets.push_back(size - j * base_offset);
        }
        base_offset *= grid_size;
    }
}
