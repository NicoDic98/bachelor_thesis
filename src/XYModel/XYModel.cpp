/**
 * @file       XYModel.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       18.05.22
 */


#include "XYModel.h"

#include <utility>

const char *XYModel::dimension_name{"dimension"};
const char *XYModel::neighbour_extent_name{"neighbour_extent"};
const char *XYModel::eta_name{"eta"};
const char *XYModel::h_name{"h"};
const char *XYModel::k_sym_name{"k_sym"};
const char *XYModel::InterpolationMatrix_name{"InterpolationMatrix"};
const char *XYModel::XYModel_name{"XYModel"};

XYModel::XYModel(double beta_, double eta_, MultiVectorX h_, int dimension_, int neighbour_extent_, int grid_size_)
        : BaseModel<MultiVectorX>(beta_, XYModel_name), eta{eta_}, h{std::move(h_)}, dimension{dimension_},
          neighbour_extent{neighbour_extent_},
          k_sym(int_pow(grid_size_, dimension_), int_pow(grid_size_, dimension_)),
          InterpolationMatrix{},
          RootModel{*this} {
    fill_connectivity_matrix(grid_size_, dimension, neighbour_extent_, k_sym);
    assert(check_internal_dimensions());
}

XYModel::XYModel(const XYModel &NewModel, InterpolationType InterpolationType_)
        : BaseModel<MultiVectorX>(NewModel), eta{NewModel.eta}, h{NewModel.h}, dimension{NewModel.dimension},
          neighbour_extent{NewModel.neighbour_extent},
          InterpolationMatrix{},
          RootModel{NewModel} {
    std::cout << "XYModel Interpolation constructor called" << std::endl;
    assert(RootModel.check_internal_dimensions());
    fill_interpolation_matrix(InterpolationType_, h[0].rows(), dimension, InterpolationMatrix);
    k_sym = InterpolationMatrix.transpose() * RootModel.k_sym * InterpolationMatrix;
    for (auto &item: h) {
        item.resize(InterpolationMatrix.cols());
    }


    //print_dimensions();
    //print_interpolation_matrix();
    assert(check_internal_dimensions());
}

XYModel::XYModel(const XYModel &NewModel)
        : BaseModel<MultiVectorX>(NewModel), eta{NewModel.eta}, h{NewModel.h}, dimension{NewModel.dimension},
          neighbour_extent{NewModel.neighbour_extent}, k_sym{NewModel.k_sym},
          InterpolationMatrix{},
          RootModel{NewModel} {
    std::cout << "XYModel copy constructor called" << std::endl;
    assert(check_internal_dimensions());
}


[[maybe_unused]] XYModel::XYModel(HighFive::Group &root)
        : BaseModel<MultiVectorX>(root, XYModel_name), eta{},
          dimension{}, neighbour_extent{}, RootModel{*this} {
    assert(name == XYModel_name);
    root.getAttribute(dimension_name).read(dimension);
    root.getAttribute(neighbour_extent_name).read(neighbour_extent);
    root.getAttribute(eta_name).read(eta);
    ReadMultiVectorX(h, root, h_name);
    ReadMatrixX(k_sym, root, k_sym_name);

    if (root.hasAttribute(InterpolationMatrix_name)) {
        ReadMatrixX(InterpolationMatrix, root, InterpolationMatrix_name);
    } else {

    }
}

void XYModel::update_pi(MultiVectorX &phi, MultiVectorX &pi, double step_size) {
    auto temp = get_force(phi);
    for (int l = 0; l < pi.size(); ++l) {
        pi[l] += step_size * temp[l];
    }
}

void XYModel::update_phi(MultiVectorX &phi, MultiVectorX &pi, double step_size) {
    for (int l = 0; l < pi.size(); ++l) {
        phi[l] += step_size * pi[l];
    }
}

double XYModel::get_action(const MultiVectorX &phi) {
    double ret{0.};
    for (int i = 0; i < phi.size(); i++) {
        ret -= eta * phi[i].dot(phi[i]);
        ret -= 0.5 * phi[i].transpose() * k_sym * phi[i];
        ret -= h[i].dot(phi[i]);
    }
    ret *= get_beta();
    return ret;
}

bool XYModel::check_dimensions(const MultiVectorX &phi) const {
    unsigned long size_a = h.size();

    long size_b = h[0].rows();
    for (const auto &elem: h) {
        if (elem.rows() != size_b) {
            return false;
        }
    }

    if (phi.size() != size_a) {
        return false;
    }
    for (const auto &elem: phi) {
        if (elem.rows() != size_b) {
            return false;
        }
    }

    if (k_sym.rows() != size_b) {
        return false;
    }
    if (k_sym.cols() != size_b) {
        return false;
    }
    return true;
}

XYModel *XYModel::get_coarser_model(InterpolationType InterpolationType_) {
    return new XYModel(*this, InterpolationType_);
}

XYModel *XYModel::get_copy_of_model() {
    return new XYModel(*this);
}

MultiVectorX XYModel::get_force(const MultiVectorX &phi) {
    auto ret = get_empty_field();
    for (int i = 0; i < phi.size(); i++) {
        ret[i] = k_sym * phi[i] + h[i];
        ret[i] += 2 * eta * phi[i];
        ret[i] *= get_beta();
    }
    return ret;
}

XYModel *XYModel::get_model_at(HighFive::Group &root) {
    return new XYModel(root);
}

void XYModel::update_fields(const MultiVectorX &phi) {
    for (int i = 0; i < h.size(); i++) {
        h[i] = (RootModel.h[i].transpose() + phi[i].transpose() * RootModel.k_sym) * InterpolationMatrix;
    }
    //TODO eta
}

void XYModel::interpolate(const MultiVectorX &phi2a, MultiVectorX &phia) {
    for (int i = 0; i < phia.size(); i++) {
        phia[i] += InterpolationMatrix * phi2a[i];
    }
}

void XYModel::pull_attributes_from_root() {
    set_beta(RootModel.get_beta());
    //TODO maybe add some more attributes to be pulled
}

MultiVectorX XYModel::get_empty_field() {
    VectorX temp(h[0].rows());
    temp.setZero();
    MultiVectorX ret;
    for (int l = 0; l < h.size(); ++l) {
        ret.push_back(temp);
    }
    return ret;
}

void XYModel::dumpToH5(HighFive::Group &root) {
    BaseModel::dumpToH5(root);
    write_static_size(dimension, root, dimension_name);
    write_static_size(neighbour_extent, root, neighbour_extent_name);
    write_static_size(eta, root, eta_name);

    WriteMultiVectorX(h, root, h_name);

    WriteMatrixX(k_sym, root, k_sym_name);

    if (InterpolationMatrix.size() != 0) {
        WriteMatrixX(InterpolationMatrix, root, InterpolationMatrix_name);
    }
}

void XYModel::load_ensemble(std::vector<MultiVectorX> &target, HighFive::DataSet &root) {
    const std::vector<size_t> shape{root.getDimensions()};
    target.resize(shape[0]);
    std::vector<std::vector<std::vector<double>>> temp;
    root.read(temp);

    for (int i = 0; i < shape[0]; ++i) {
        for (auto &item: temp[i]) {
            target[i].push_back(VectorX::Map(&item[0], item.size()));
        }
    }
}

HighFive::DataSet
XYModel::dump_ensemble(std::vector<MultiVectorX> &target, HighFive::Group &root, std::string sub_name) {
    std::vector<size_t> offset;
    HighFive::DataSet target_dataset = add_to_expandable_dataset(
            root, sub_name,
            {target.size(), target[0].size(), static_cast<unsigned long>(target[0][0].rows())},
            offset);

    std::vector<size_t> count{1, 1, static_cast<unsigned long>(target[0][0].rows())};

    auto offset1 = offset[1];
    for (auto &elem: target) {
        offset[1] = offset1;
        for (auto &inner_elem: elem) {
            target_dataset.select(offset, count).write_raw(inner_elem.data());
            offset[1]++;
        }
        offset[0]++;
    }

    return target_dataset;
}

double XYModel::get_artificial_energy(const MultiVectorX &phi, const MultiVectorX &pi) {
    double ret{0.};
    for (const auto &elem: pi) {
        ret += elem.dot(elem);
    }
    ret *= 0.5;
    return ret + get_action(phi);
}

MultiVectorX XYModel::get_pi(std::default_random_engine &generator) {
    auto ret = get_empty_field();
    std::normal_distribution<double> gauss(0, 1);
    for (auto &elem: ret) {
        for (auto &inner_elem: elem) {
            inner_elem = gauss(generator);
        }
    }
    return ret;
}

void XYModel::renormalize(MultiVectorX &phi) {
    for (int l = 0; l < phi[0].rows(); ++l) {
        double norm{0.};
        for (auto elem: phi) {
            norm += elem[l] * elem[l];
        }
        norm = sqrt(norm);
        for (auto &elem: phi) {
            elem[l] /= norm;
        }
    }
}

double XYModel::get_vector_length_squared(const MultiVectorX &phi) {
    double norm{0.};
    for (const auto &elem: phi) {
        norm += elem.dot(elem);
    }
    return norm / static_cast<double>(phi[0].rows());
}




