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
        : BaseModel<VectorX>(beta_), sqrt_beta{sqrt(beta_)}, h{std::move(h_)}, eta{std::move(eta_)},
          dimension{dimension_},
          k_sym(int_pow(grid_size_, dimension_), int_pow(grid_size_, dimension_)),
          k_rec(int_pow(grid_size_, dimension_), int_pow(grid_size_, dimension_)),
          InterpolationMatrix{} {

    fill_connectivity_matrix(neighbour_extent_, grid_size_);
    add_offset_to_connectivity_matrix(offset_);

    k_rec = k_sym;
    assert(check_internal_dimensions());
}

IsingModel::IsingModel(const IsingModel &NewModel, InterpolationType InterpolationType_)
        : BaseModel<VectorX>(NewModel.get_beta()), sqrt_beta{sqrt(NewModel.get_beta())}, h{NewModel.h},
          eta{NewModel.eta}, dimension{NewModel.dimension}, InterpolationMatrix{} {
    std::cout << "Hello" << std::endl;
    assert(NewModel.check_internal_dimensions());
    int resize_factor = fill_interpolation_matrix(InterpolationType_);
    k_sym = InterpolationMatrix.transpose() * NewModel.k_sym * InterpolationMatrix;
    k_rec = NewModel.k_rec * InterpolationMatrix;
    h.resize(h.rows() / resize_factor);

    print_dimensions();
    assert(check_internal_dimensions());
}

IsingModel::IsingModel(const IsingModel &NewModel)
        : BaseModel<VectorX>(NewModel.get_beta()), sqrt_beta{sqrt(NewModel.get_beta())}, h{NewModel.h},
          eta{NewModel.eta}, dimension{NewModel.dimension}, k_sym{NewModel.k_sym}, k_rec{NewModel.k_rec},
          InterpolationMatrix{NewModel.InterpolationMatrix} {
    std::cout << "hi" << std::endl;

    print_dimensions();
    assert(check_internal_dimensions());
}

double IsingModel::get_action(const VectorX &phi) {
    VectorX var_phi{k_rec * phi};
    double ret{0.};
    for (int i = 0; i < eta.rows(); ++i) {
        ret -= log(cosh(eta(i) + sqrt_beta * var_phi(i)) / cosh(eta(i)));
    }
    ret += 0.5 * phi.transpose() * k_sym * phi;
    ret -= sqrt_beta * h.dot(phi);
    return ret;
}

VectorX IsingModel::get_force(const VectorX &phi) {
    VectorX var_phi1{k_rec * phi};
    VectorX var_phi2{k_sym * phi};
    VectorX temp{eta + sqrt_beta * var_phi1};
    for (auto &elem: temp) {
        elem = tanh(elem);
    }
    return -var_phi2 + sqrt_beta * h + sqrt_beta * k_rec.transpose() * temp;
}

double IsingModel::get_magnetization(const VectorX &phi) {

    return (phi / sqrt_beta - (k_sym.inverse()) * h).sum() / phi.rows();
}

void IsingModel::fill_connectivity_matrix(int neighbour_extent, int grid_size) {
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

void IsingModel::add_offset_to_connectivity_matrix(double offset) {
    MatrixX OffsetMatrix(k_sym.rows(), k_sym.cols());
    OffsetMatrix.setIdentity();
    k_sym += offset * OffsetMatrix;
}

void IsingModel::print_connectivity_matrix() {
    std::cout << k_sym << std::endl;
}

bool IsingModel::check_dimensions(const VectorX &phi) const {
    long size_a = h.rows();
    long size_base = eta.rows();

    if (phi.rows() != size_a) {
        return false;
    }

    if (k_sym.rows() != size_a) {
        return false;
    }
    if (k_sym.cols() != size_a) {
        return false;
    }

    if (k_rec.rows() != size_base) {
        return false;
    }
    if (k_rec.cols() != size_a) {
        return false;
    }
    return true;
}

void IsingModel::print_name() {
    std::cout << "I'm the IsingModel" << std::endl;
}

IsingModel *IsingModel::get_coarser_model(InterpolationType InterpolationType_) {
    return new IsingModel(*this, InterpolationType_);
}

IsingModel *IsingModel::get_copy_of_model() {
    return new IsingModel(*this);
}

int IsingModel::fill_interpolation_matrix(InterpolationType InterpolationType_) {
    int resize_factor{1};
    switch (InterpolationType_) {
        case InterpolationType::Checkerboard:
            resize_factor = 4;
            InterpolationMatrix.resize(h.rows(), h.rows() / resize_factor);//TODO What about 6/4=1?
            break;
        case InterpolationType::Black_White:
            resize_factor = 2;
            InterpolationMatrix.resize(h.rows(), h.rows() / resize_factor);
            break;
    }
    return resize_factor;
}

void IsingModel::print_dimensions() {
    std::cout << "h dimensions:\t (" << h.rows() << ", " << h.cols() << ')' << std::endl;
    std::cout << "eta dimensions:\t (" << eta.rows() << ", " << eta.cols() << ')' << std::endl;
    std::cout << "k_sym dimensions:\t (" << k_sym.rows() << ", " << k_sym.cols() << ')' << std::endl;
    std::cout << "k_rec dimensions:\t (" << k_rec.rows() << ", " << k_rec.cols() << ')' << std::endl;
}






