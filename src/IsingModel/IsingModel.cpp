/**
 * @file       IsingModel.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       23.03.22
 */


#include <iostream>
#include "IsingModel.h"

inline int int_pow(int x, int y) {
    assert(y > 0);
    int result = 1;
    for (int i = 0; i < y; ++i) {
        result *= x;
    }
    return result;
}

IsingModel::IsingModel(int dimension_, int grid_size_, int neighbour_extent_, double offset)
        : dimension{dimension_}, grid_size{grid_size_}, neighbour_extent{neighbour_extent_},
          connectivity_matrix(int_pow(grid_size_, dimension_), int_pow(grid_size_, dimension_)) {
    fill_connectivity_matrix();
    add_offset_to_connectivity_matrix(offset);
}

double IsingModel::get_action(const VectorX &phi) {
    return 0;
}

VectorX IsingModel::get_force(const VectorX &phi) {
    return VectorX();
}

double IsingModel::get_magnetization(const VectorX &phi) {
    return 0;
}

void IsingModel::fill_connectivity_matrix() {
    connectivity_matrix.setZero();

    std::vector<long> index_offsets;
    long size = connectivity_matrix.rows();
    generate_index_offsets(index_offsets, size);

    for (auto elem: index_offsets) {
        std::cout << elem << '\n';
    }

    for (long m = 0; m < size; ++m) {
        for (auto index_offset: index_offsets) {
            connectivity_matrix(m, (m + index_offset) % size) += 1;
        }

    }
}

void IsingModel::add_offset_to_connectivity_matrix(double offset) {
    MatrixX OffsetMatrix(connectivity_matrix.rows(), connectivity_matrix.cols());
    OffsetMatrix.setIdentity();
    connectivity_matrix += offset * OffsetMatrix;
}

void IsingModel::print_connectivity_matrix() {
    std::cout << connectivity_matrix;
}

void IsingModel::generate_index_offsets(std::vector<long> &index_offsets, long size) {
    long base_offset{1};
    for (int i = 0; i < dimension; ++i) {
        for (long j = 1; j <= neighbour_extent; ++j) {
            index_offsets.push_back(j * base_offset);
            index_offsets.push_back(size - j * base_offset);
        }
        base_offset *= grid_size;
    }
}
