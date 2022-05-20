/**
 * @file       XYModel.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       18.05.22
 */


#include "XYModel.h"

const char *XYModel::dimension_name{"dimension"};
const char *XYModel::neighbour_extent_name{"neighbour_extent"};
const char *XYModel::h_name{"h"};
const char *XYModel::k_sym_name{"k_sym"};
const char *XYModel::InterpolationMatrix_name{"InterpolationMatrix"};
const char *XYModel::XYModel_name{"XYModel"};

XYModel::XYModel(double beta_, const MultiVectorX &h_, int dimension_, int neighbour_extent_, int grid_size_)
        : BaseModel<MultiVectorX>(beta_, XYModel_name), h{h_}, dimension{dimension_},
          neighbour_extent{neighbour_extent_},
          k_sym(int_pow(grid_size_, dimension_), int_pow(grid_size_, dimension_)),
          InterpolationMatrix{},
          RootModel{*this} {
    fill_connectivity_matrix(grid_size_, dimension, neighbour_extent_, k_sym);
    assert(check_internal_dimensions());
}

XYModel::XYModel(const XYModel &NewModel, InterpolationType InterpolationType_)
        : BaseModel<MultiVectorX>(NewModel), h{NewModel.h}, dimension{NewModel.dimension},
          neighbour_extent{NewModel.neighbour_extent},
          InterpolationMatrix{},
          RootModel{NewModel} {
    std::cout << "XYModel Interpolation constructor called" << std::endl;
    assert(RootModel.check_internal_dimensions());
    fill_interpolation_matrix(InterpolationType_, h[0].rows(), dimension, InterpolationMatrix);
    k_sym = InterpolationMatrix.transpose() * RootModel.k_sym * InterpolationMatrix;
    h.resize(InterpolationMatrix.cols());

    //print_dimensions();
    //print_interpolation_matrix();
    assert(check_internal_dimensions());
}

XYModel::XYModel(const XYModel &NewModel)
        : BaseModel<MultiVectorX>(NewModel), h{NewModel.h}, dimension{NewModel.dimension},
          neighbour_extent{NewModel.neighbour_extent}, k_sym{NewModel.k_sym},
          InterpolationMatrix{},
          RootModel{NewModel} {
    std::cout << "XYModel copy constructor called" << std::endl;
    assert(check_internal_dimensions());
}


[[maybe_unused]] XYModel::XYModel(HighFive::Group &root)
        : BaseModel<MultiVectorX>(root, XYModel_name),
          dimension{}, neighbour_extent{}, RootModel{*this} {
    assert(name == XYModel_name);
    root.getAttribute(dimension_name).read(dimension);
    root.getAttribute(neighbour_extent_name).read(neighbour_extent);
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
        ret -= 0.5 * phi[i].transpose() * k_sym * phi[i];
        ret -= h[i].dot(phi[i]);
    }
    ret *= get_beta();
    return ret;
}

bool XYModel::check_dimensions(const MultiVectorX &phi) const {
    return false;
}

XYModel *XYModel::get_coarser_model(InterpolationType InterpolationType_) {
    return nullptr;
}

XYModel *XYModel::get_copy_of_model() {
    return nullptr;
}

MultiVectorX XYModel::get_force(const MultiVectorX &phi) {
    return MultiVectorX();
}

XYModel *XYModel::get_model_at(HighFive::Group &root) {
    return nullptr;
}

void XYModel::update_fields(const MultiVectorX &phi) {

}

void XYModel::interpolate(const MultiVectorX &phi2a, MultiVectorX &phia) {

}

void XYModel::pull_attributes_from_root() {

}

MultiVectorX XYModel::get_empty_field() {
    return MultiVectorX();
}

void XYModel::dumpToH5(HighFive::Group &root) {
    BaseModel::dumpToH5(root);
}

void XYModel::load_ensemble(std::vector<MultiVectorX> &target, HighFive::DataSet &root) {

}

HighFive::DataSet
XYModel::dump_ensemble(std::vector<MultiVectorX> &target, HighFive::Group &root, std::string sub_name) {
    return HighFive::DataSet();
}




