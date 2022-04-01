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
          dimension{dimension_}, grid_side_length{grid_size_},
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
          eta{NewModel.eta}, dimension{NewModel.dimension}, grid_side_length{NewModel.grid_side_length},
          InterpolationMatrix{} {
    std::cout << "Hello" << std::endl;
    assert(NewModel.check_internal_dimensions());
    grid_side_length = fill_interpolation_matrix(InterpolationType_, h.rows(), grid_side_length);
    k_sym = InterpolationMatrix.transpose() * NewModel.k_sym * InterpolationMatrix;
    k_rec = NewModel.k_rec * InterpolationMatrix;
    h.resize(int_pow(grid_side_length, dimension));

    print_dimensions();
    print_interpolation_matrix();
    assert(check_internal_dimensions());
}

IsingModel::IsingModel(const IsingModel &NewModel)
        : BaseModel<VectorX>(NewModel.get_beta()), sqrt_beta{sqrt(NewModel.get_beta())}, h{NewModel.h},
          eta{NewModel.eta}, dimension{NewModel.dimension}, grid_side_length{NewModel.grid_side_length},
          k_sym{NewModel.k_sym}, k_rec{NewModel.k_rec},
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

int
IsingModel::fill_interpolation_matrix(InterpolationType InterpolationType_, long fine_size, int fine_grid_side_length) {
    int coarse_grid_side_length{1};
    int coarse_lambda{1};
    switch (InterpolationType_) {
        case InterpolationType::Checkerboard:

            if (fine_grid_side_length % 2 == 0) {
                coarse_grid_side_length = fine_grid_side_length / 2;
            } else {
                coarse_grid_side_length = (fine_grid_side_length + 1) / 2;
                //TODO: think more about this case
            }

            coarse_lambda = int_pow(coarse_grid_side_length, dimension);

            InterpolationMatrix.resize(fine_size, coarse_lambda);
            InterpolationMatrix.setZero();


            for (long m = 0; m < InterpolationMatrix.rows(); ++m) {
                long base_offset{1};
                long coarse_offset{1};
                double factor{1};
                std::vector<long> fine_interpolation_ids{m};
                for (int i = 0; i < dimension; ++i) {
                    if (((m / base_offset) % fine_grid_side_length) % 2 == 1) {
                        factor *= 0.5;
                        std::vector<long> temp{fine_interpolation_ids};
                        fine_interpolation_ids.clear();
                        for (long elem: temp) {
                            fine_interpolation_ids.push_back(
                                    (elem + base_offset) % (base_offset * fine_grid_side_length) +
                                    (base_offset * fine_grid_side_length) *
                                    (elem / (base_offset * fine_grid_side_length)));
                            fine_interpolation_ids.push_back(
                                    (elem - base_offset + (base_offset * fine_grid_side_length)) %
                                    (base_offset * fine_grid_side_length) +
                                    (base_offset * fine_grid_side_length) *
                                    (elem / (base_offset * fine_grid_side_length)));
                        }
                    }

                    base_offset *= fine_grid_side_length;
                }
                //assert(factor == fine_interpolation_ids.size());
                for (auto elem: fine_interpolation_ids) {
                    base_offset = 1;
                    coarse_offset = 1;
                    long coarse_id{0};
                    for (int i = 0; i < dimension; ++i) {
                        coarse_id +=
                                ((elem % (base_offset * fine_grid_side_length)) / (base_offset * 2)) * coarse_offset;

                        base_offset *= fine_grid_side_length;
                        coarse_offset *= coarse_grid_side_length;
                    }
                    InterpolationMatrix(m, coarse_id) = factor;
                }

            }
            break;
        case InterpolationType::Black_White:
            coarse_grid_side_length = 42;//TODO think if this can even be generalized to e.g. 3 dimensions
            InterpolationMatrix.resize(fine_size, coarse_grid_side_length);
            break;
    }
    return coarse_grid_side_length;
}

void IsingModel::print_dimensions() {
    std::cout << "h dimensions:\t (" << h.rows() << ", " << h.cols() << ')' << std::endl;
    std::cout << "eta dimensions:\t (" << eta.rows() << ", " << eta.cols() << ')' << std::endl;
    std::cout << "k_sym dimensions:\t (" << k_sym.rows() << ", " << k_sym.cols() << ')' << std::endl;
    std::cout << "k_rec dimensions:\t (" << k_rec.rows() << ", " << k_rec.cols() << ')' << std::endl;
    std::cout << "InterpolationMatrix dimensions:\t (" << InterpolationMatrix.rows() << ", " << InterpolationMatrix.cols() << ')' << std::endl;
}

void IsingModel::print_interpolation_matrix() {
    std::cout << InterpolationMatrix << std::endl;
}






