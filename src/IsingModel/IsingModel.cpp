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
    VectorX var_phi{k_rec * phi};
    double ret{0.};
    std::cout << '\n' << phi << '\n' << std::endl;
    for (int i = 0; i < eta.rows(); ++i) {
        ret -= log(cosh(eta(i) + sqrt_beta * var_phi(i)) / cosh(eta(i)));
    }
    ret += 0.5 * phi.transpose() * k_sym * phi;
    ret -= sqrt_beta * h.dot(phi);
    return ret;
}

VectorX IsingModel::get_force(const VectorX &phi) {
    VectorX var_phi{k_rec * phi};
    VectorX temp{eta + sqrt_beta * var_phi};
    for (auto &elem: temp) {
        elem = tanh(elem);
    }
    //std::cout << "2:\n" << -var_phi + sqrt_beta * h + sqrt_beta * k_rec.transpose() * temp << std::endl;
    return -var_phi + sqrt_beta * h + sqrt_beta * k_rec.transpose() * temp;
}

double IsingModel::get_magnetization(const VectorX &phi) {
    return 0;
}

void IsingModel::fill_connectivity_matrix(int dimension, int neighbour_extent, int grid_size) {
    k_sym.setZero();

    long lambda = k_sym.rows();

    for (long m = 0; m < lambda; ++m) {
        long base_offset{1};
        for (int i = 0; i < dimension; ++i) {
            for (long j = 1; j <= neighbour_extent; ++j) {
                k_sym(m, (m + base_offset * neighbour_extent) % (base_offset * grid_size)
                         + (base_offset * grid_size) * (m / (base_offset * grid_size))) += 1;
                // this additional + is needed in cas of eg ((0,1),(2,3)) to get the correct connections for 2
                k_sym(m, (m - base_offset * neighbour_extent + (base_offset * grid_size)) % (base_offset * grid_size)
                         + (base_offset * grid_size) * (m / (base_offset * grid_size))) += 1;
            }
            base_offset *= grid_size;
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
