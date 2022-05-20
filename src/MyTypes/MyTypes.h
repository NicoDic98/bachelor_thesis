/**
 * @file       MyTypes.h
 * @brief      Definition of used types
 * @author     nico
 * @version    0.0.2
 * @date       24.03.22
 */


#ifndef BACHELOR_THESIS_MYTYPES_H
#define BACHELOR_THESIS_MYTYPES_H

#include <eigen3/Eigen/Dense>
#include <highfive/H5File.hpp>

/**
  * @brief Vector type definition
  */
typedef Eigen::VectorXd VectorX;

/**
 * @brief Matrix type definition
 */
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixX;

typedef std::vector<VectorX> MultiVectorX;

/**
 * @brief Interpolation types
 */
enum class InterpolationType {
    Checkerboard, Black_White
};


/**
 * @brief Simple Power function to calculate the powers of integers
 * @param x Base
 * @param y Exponent
 * @return Base^Exponent
 */
inline int int_pow(int x, int y) {
    assert(y >= 0);
    int result = 1;
    for (int i = 0; i < y; ++i) {
        result *= x;
    }
    return result;
}

/**
 * @brief Fills the connectivity matrix \a k_sym for the given hyper cube of dimension \a dimension and side length \p grid_size
 * @param grid_size Side length
 */
void fill_connectivity_matrix(int grid_size, int dimension, int neighbour_extent, MatrixX &k_sym);

/**
 * @brief Fills the \a InterpolationMatrix
 * @param InterpolationType_ Type of interpolation to be used
 * @param fine_size Finer grid total size
 */
void fill_interpolation_matrix(InterpolationType InterpolationType_, long fine_size, int dimension,
                               MatrixX &InterpolationMatrix);

/**
 * @brief Writes \p vector_to_write to \p root as dataset with the name \p name
 * @param vector_to_write Vector to be written
 * @param root Group to write to
 * @param name Name of the target dataset
 */
void WriteVectorX(VectorX &vector_to_write, HighFive::Group &root, std::string name);

/**
 * @brief Reads \p vector_to_read_into from \p root from the dataset with the name \p name
 * @param vector_to_read_into Vector to read into
 * @param root Group to read from
 * @param name Name of the target dataset
 */
void ReadVectorX(VectorX &vector_to_read_into, HighFive::Group &root, std::string name);

/**
 * @brief Writes \p matrix_to_write to \p root as dataset with the name \p name
 * @param matrix_to_write Matrix to be written
 * @param root Group to write to
 * @param name Name of the target dataset
 */
void WriteMatrixX(MatrixX &matrix_to_write, HighFive::Group &root, std::string name);

/**
 * @brief Reads \p matrix_to_read_into from \p root from the dataset with the name \p name
 * @param matrix_to_read_into Matrix to read into
 * @param root Group to read from
 * @param name Name of the target dataset
 */
void ReadMatrixX(MatrixX &matrix_to_read_into, HighFive::Group &root, std::string name);

/**
 * @brief Writes constant size \p var to \p root as attribute with the name \p name
 * @param var Variable to be written
 * @param root Group/DataSet to write to
 * @param name Name of the target Attribute
 */
template<typename T, typename group_dataset>
void write_static_size(T var, HighFive::AnnotateTraits<group_dataset> &root, std::string name) {
    if (root.hasAttribute(name)) {
        root.getAttribute(name).write(var);
        //type should always be double, no size checks
    } else {
        root.createAttribute(name, var);
    }
}

/**
 * @brief Adds/Creates the dataset \a sub_name under \a root with dimensions \a dims
 * @param root Group under which the dataset lies
 * @param sub_name Name of the dataset
 * @param dims Dimensions of the dataset to be added
 * @param offset Offset to be used in the following writing task
 * @param override Rather to overwrite the existing dataset or not
 * @param dims_to_resize Dimensions which should be/are resizeable
 * @return Dataset into which a write-task can be done using \a offset
 */
HighFive::DataSet add_to_expandable_dataset(HighFive::Group &root, const std::string &sub_name,
                                            const std::vector<size_t> &dims, std::vector<size_t> &offset,
                                            bool override = false, const std::vector<size_t> &dims_to_resize = {0});

#endif //BACHELOR_THESIS_MYTYPES_H
