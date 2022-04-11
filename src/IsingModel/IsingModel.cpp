/**
 * @file       IsingModel.cpp
 * @brief      Definition of the Ising Model methods
 * @author     nico
 * @version    0.0.1
 * @date       23.03.22
 */


#include <iostream>
#include <utility>
#include "IsingModel.h"


const char *IsingModel::dimension_name{"dimension"};
const char *IsingModel::grid_side_length_name{"grid_side_length"};
const char *IsingModel::neighbour_extent_name{"neighbour_extent"};
const char *IsingModel::h_name{"h"};
const char *IsingModel::eta_name{"eta"};
const char *IsingModel::connectivity_offset_name{"connectivity_offset"};
const char *IsingModel::k_sym_name{"k_sym"};
const char *IsingModel::k_rec_name{"k_rec"};
const char *IsingModel::InterpolationMatrix_name{"InterpolationMatrix"};
const char *IsingModel::IsingModel_name{"IsingModel"};


IsingModel::IsingModel(double beta_, VectorX h_, VectorX eta_, double offset_, int dimension_,
                       int neighbour_extent_, int grid_size_)
        : BaseModel<VectorX>(beta_, IsingModel_name), sqrt_beta{sqrt(beta_)}, h{std::move(h_)},
          eta{std::move(eta_)}, connectivity_offset{offset_},
          dimension{dimension_}, neighbour_extent{neighbour_extent_}, grid_side_length{grid_size_},
          k_sym(int_pow(grid_size_, dimension_), int_pow(grid_size_, dimension_)),
          k_rec(int_pow(grid_size_, dimension_), int_pow(grid_size_, dimension_)),
          InterpolationMatrix{}, FinerModel{*this} {

    fill_connectivity_matrix(neighbour_extent_, grid_size_);
    add_offset_to_connectivity_matrix(offset_ + (2 * neighbour_extent_ * dimension));
    //todo: test that this really yields a positive definit connectivity matrix (search for min eigenvalue before adding the offset)
    k_rec = k_sym;
    assert(check_internal_dimensions());
}

IsingModel::IsingModel(const IsingModel &NewModel, InterpolationType InterpolationType_)
        : BaseModel<VectorX>(NewModel), sqrt_beta{sqrt(NewModel.get_beta())}, h{NewModel.h},
          eta{NewModel.eta}, connectivity_offset{NewModel.connectivity_offset}, dimension{NewModel.dimension},
          neighbour_extent{NewModel.neighbour_extent}, grid_side_length{NewModel.grid_side_length},
          InterpolationMatrix{}, FinerModel{NewModel} {
    std::cout << "Isingmodel Interpolation constructor called" << std::endl;
    assert(FinerModel.check_internal_dimensions());
    grid_side_length = fill_interpolation_matrix(InterpolationType_, h.rows(), grid_side_length);
    set_k_sym(InterpolationMatrix.transpose() * FinerModel.k_sym * InterpolationMatrix);
    k_rec = FinerModel.k_rec * InterpolationMatrix;
    h.resize(int_pow(grid_side_length, dimension));

    //print_dimensions();
    //print_interpolation_matrix();
    assert(check_internal_dimensions());
}

IsingModel::IsingModel(const IsingModel &NewModel)
        : BaseModel<VectorX>(NewModel), sqrt_beta{sqrt(NewModel.get_beta())}, h{NewModel.h},
          eta{NewModel.eta}, connectivity_offset{NewModel.connectivity_offset}, dimension{NewModel.dimension},
          neighbour_extent{NewModel.neighbour_extent}, grid_side_length{NewModel.grid_side_length},
          k_sym(NewModel.k_sym.rows(), NewModel.k_sym.cols()), k_rec{NewModel.k_rec},
          InterpolationMatrix{NewModel.InterpolationMatrix}, FinerModel{*this} {
    std::cout << "Isingmodel copy constructor called" << std::endl;
    set_k_sym(NewModel.k_sym);
    //print_dimensions();
    assert(check_internal_dimensions());
}

[[maybe_unused]] IsingModel::IsingModel(HighFive::Group &root)
        : BaseModel<VectorX>(root, IsingModel_name),
          sqrt_beta{sqrt(get_beta())}, dimension{}, neighbour_extent{}, grid_side_length{}, connectivity_offset{},
          FinerModel{*this} {
    assert(name == IsingModel_name);
    root.getAttribute(dimension_name).read(dimension);
    root.getAttribute(neighbour_extent_name).read(neighbour_extent);
    root.getAttribute(grid_side_length_name).read(grid_side_length);
    ReadVectorX(h, root, h_name);
    ReadVectorX(eta, root, eta_name);
    root.getAttribute(connectivity_offset_name).read(connectivity_offset);
    ReadMatrixX(k_sym, root, k_sym_name);
    set_k_sym(k_sym);
    ReadMatrixX(k_rec, root, k_rec_name);

    if (root.hasAttribute(InterpolationMatrix_name)) {
        ReadMatrixX(InterpolationMatrix, root, InterpolationMatrix_name);
    } else {

    }
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

    return (phi / sqrt_beta - (k_sym_inverse) * h).sum() / static_cast<double>(phi.rows());
}

double IsingModel::get_magnetization_squared(const VectorX &phi) {

    double m_squared{0.};
    VectorX k_sym_inv_h{k_sym_inverse * h};
    MatrixX phi_phi(phi.rows(), phi.rows());
    MatrixX phi_k_sym_inv_h(phi.rows(), phi.rows());
    MatrixX k_sym_inv_h_k_sym_inv_h(phi.rows(), phi.rows());


    for (int i = 0; i < phi.rows(); ++i) {
        for (int j = 0; j < phi.rows(); ++j) {
            phi_phi(i, j) = phi(i) * phi(j);
            phi_k_sym_inv_h(i, j) = phi(i) * k_sym_inv_h(j);
            k_sym_inv_h_k_sym_inv_h(i, j) = k_sym_inv_h(i) * k_sym_inv_h(j);
        }
    }
    m_squared -= k_sym_inverse.sum() / get_beta();
    m_squared += phi_phi.sum() / get_beta();
    m_squared -= 2 * phi_k_sym_inv_h.sum() / sqrt_beta;
    m_squared += k_sym_inv_h_k_sym_inv_h.sum();

    return m_squared / static_cast<double>(phi.rows() * phi.rows());
}

double IsingModel::get_energy(const VectorX &phi) {
    double e{0.};

    e += 0.5 * h.transpose() * k_sym_inverse * h;

    e -= 0.5 / sqrt_beta * h.transpose() * phi;

    VectorX var_phi{k_rec * phi};
    VectorX temp{eta + sqrt_beta * var_phi};
    for (int i = 0; i < var_phi.rows(); ++i) {
        temp(i) = tanh(temp(i)) * var_phi(i);
    }
    e -= 0.5 * temp.sum() / sqrt_beta;

    e /= static_cast<double >(phi.rows());
    e += 0.5 * connectivity_offset;
    return e;
}

double IsingModel::get_energy_squared(const VectorX &phi) {
    double e_squared{0.};

    e_squared += 1. / (sqrt_beta) * h.transpose() * phi;

    VectorX var_phi{k_rec * phi};
    VectorX temp{eta + sqrt_beta * var_phi};
    for (int i = 0; i < var_phi.rows(); ++i) {
        temp(i) = tanh(temp(i)) * var_phi(i) / sqrt_beta
                  - var_phi(i) * var_phi(i) / pow(cos(temp(i)), 2);
    }
    e_squared += temp.sum();

    e_squared *= 0.25 / (static_cast<double>(phi.rows()) * get_beta());

    e_squared += pow(get_energy(phi), 2);

    return e_squared;
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
    set_k_sym(k_sym);
}

void IsingModel::add_offset_to_connectivity_matrix(double offset) {
    MatrixX OffsetMatrix(k_sym.rows(), k_sym.cols());
    OffsetMatrix.setIdentity();
    set_k_sym(k_sym + offset * OffsetMatrix);
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
    std::cout << name << std::endl;
}

IsingModel *IsingModel::get_coarser_model(InterpolationType InterpolationType_) {
    return new IsingModel(*this, InterpolationType_);
}

IsingModel *IsingModel::get_copy_of_model() {
    return new IsingModel(*this);
}

IsingModel *IsingModel::get_model_at(HighFive::Group &root) {
    return new IsingModel(root);
}

int
IsingModel::fill_interpolation_matrix(InterpolationType InterpolationType_, long fine_size, int fine_grid_side_length) {
    int coarse_grid_side_length{1};
    int coarse_lambda;
    switch (InterpolationType_) {
        case InterpolationType::Checkerboard:

            /*
             * Compute new sizes
             */
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
                long coarse_offset;
                double factor{1};
                /*
                 * Get ids of fine level sources of the current m
                 */
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

                /*
                 * Fill interpolation matrix, by calculating the coarse level id of each source
                 */
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

[[maybe_unused]] void IsingModel::print_dimensions() {
    std::cout << h_name << " dimensions:\t (" << h.rows() << ", " << h.cols() << ')' << std::endl;
    std::cout << eta_name << " dimensions:\t (" << eta.rows() << ", " << eta.cols() << ')' << std::endl;
    std::cout << k_sym_name << " dimensions:\t (" << k_sym.rows() << ", " << k_sym.cols() << ')' << std::endl;
    std::cout << k_rec_name << " dimensions:\t (" << k_rec.rows() << ", " << k_rec.cols() << ')' << std::endl;
    std::cout << InterpolationMatrix_name << " dimensions:\t (" << InterpolationMatrix.rows() << ", "
              << InterpolationMatrix.cols() << ')' << std::endl;
}

[[maybe_unused]] void IsingModel::print_interpolation_matrix() {
    std::cout << InterpolationMatrix << std::endl;
}

void IsingModel::update_fields(const VectorX &phi) {
    h = (FinerModel.h.transpose() - phi.transpose() * FinerModel.k_sym / sqrt_beta) * InterpolationMatrix;
    eta = FinerModel.eta + sqrt_beta * FinerModel.k_rec * phi;
}

void IsingModel::interpolate(const VectorX &phi2a, VectorX &phia) {
    phia += InterpolationMatrix * phi2a;
}

VectorX IsingModel::get_empty_field() {
    VectorX temp(h.rows());
    temp.setZero();
    return temp;
}

void IsingModel::dumpToH5(HighFive::Group &root) {
    BaseModel<VectorX>::dumpToH5(root);
    if (root.hasAttribute(dimension_name)) {
        HighFive::Attribute temp = root.getAttribute(dimension_name);
        temp.write(dimension);
    } else {
        root.createAttribute(dimension_name, dimension);
    }

    if (root.hasAttribute(neighbour_extent_name)) {
        HighFive::Attribute temp = root.getAttribute(neighbour_extent_name);
        temp.write(neighbour_extent);
    } else {
        root.createAttribute(neighbour_extent_name, neighbour_extent);
    }

    if (root.hasAttribute(grid_side_length_name)) {
        HighFive::Attribute temp = root.getAttribute(grid_side_length_name);
        temp.write(grid_side_length);
    } else {
        root.createAttribute(grid_side_length_name, grid_side_length);
    }

    WriteVectorX(h, root, h_name);

    WriteVectorX(eta, root, eta_name);

    if (root.hasAttribute(connectivity_offset_name)) {
        HighFive::Attribute temp = root.getAttribute(connectivity_offset_name);
        temp.write(connectivity_offset);
    } else {
        root.createAttribute(connectivity_offset_name, connectivity_offset);
    }

    WriteMatrixX(k_sym, root, k_sym_name);
    WriteMatrixX(k_rec, root, k_rec_name);

    if (InterpolationMatrix.size() != 0) {
        WriteMatrixX(InterpolationMatrix, root, InterpolationMatrix_name);
    }

}

void IsingModel::pull_attributes_from_finer_level() {
    set_beta(FinerModel.get_beta());
    //TODO maybe add some more attributes to be pulled
}

void IsingModel::set_k_sym(const MatrixX &k_sym_new) {
    k_sym = k_sym_new;
    k_sym_inverse = k_sym.inverse();
}

void IsingModel::load_ensemble(std::vector<VectorX> &target, HighFive::DataSet &root) {
    const std::vector<size_t> shape{root.getDimensions()};
    target.resize(shape[0]);
    std::vector<std::vector<double>> temp;
    root.read(temp);

    for (int i = 0; i < shape[0]; ++i) {
        VectorX vec_temp(shape[1]);
        vec_temp = VectorX::Map(&temp[i][0], vec_temp.size());
        target[i] = vec_temp;
    }
}

HighFive::DataSet IsingModel::dump_ensemble(std::vector<VectorX> &target, HighFive::Group &root, std::string sub_name) {
    HighFive::DataSet target_dataset;
    std::vector<size_t> offset{0, 0};
    std::vector<size_t> count{1, 0};
    if (root.exist(sub_name)) {
        target_dataset = root.getDataSet(sub_name);
        const std::vector<size_t> shape{target_dataset.getDimensions()};

        if (shape.size() == 2) {
            if (shape[1] != target[0].rows()) {
                std::cerr << "Vector sizes not matching\n";
                exit(-1);
            }
            offset[0] = shape[0];
            count[1] = shape[1];
            target_dataset.resize({shape[0] + target.size(), shape[1]});

        } else if (shape.size() == 3) {
            if (shape[2] != 1) {
                std::cerr << "Shapes not matching\n";
                exit(-1);
            }
            if (shape[1] != target[0].rows()) {
                std::cerr << "Vector sizes not matching\n";
                exit(-1);
            }
            //TODO

        } else {
            std::cerr << "Can't extend dataset\n";
            exit(-1);
        }

    } else {
        auto my_dataspace = HighFive::DataSpace(
                {target.size(), static_cast<unsigned long>(target[0].rows())},
                {HighFive::DataSpace::UNLIMITED, static_cast<unsigned long>(target[0].rows())});
        HighFive::DataSetCreateProps props;
        props.add(HighFive::Chunking({1, static_cast<unsigned long>(target[0].rows())}));
        target_dataset = root.createDataSet(sub_name, my_dataspace,
                                            HighFive::create_datatype<double>(), props);
        count[1] = target[0].rows();
    }
    for (auto & elem : target) {
        target_dataset.select(offset,count).write_raw(elem.data());
        offset[0]++;
    }

    return target_dataset;
}
