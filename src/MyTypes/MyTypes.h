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

/**
 * @brief Interpolation types
 */
enum class InterpolationType {
    Checkerboard, Black_White
};

/**
 * @brief Writes \p vector_to_write to \p root as attribute with the name \p name
 * @param vector_to_write Vector to be written
 * @param root Group to write to
 * @param name Name of the target Attribute
 */
void WriteVectorX(VectorX &vector_to_write, HighFive::Group &root, std::string name);

/**
 * @brief Reads \p vector_to_read_into from \p root from the attribute with the name \p name
 * @param vector_to_read_into Vector to read into
 * @param root Group to read from
 * @param name Name of the target Attribute
 */
void ReadVectorX(VectorX &vector_to_read_into, HighFive::Group &root, std::string name);

/**
 * @brief Writes \p matrix_to_write to \p root as attribute with the name \p name
 * @param matrix_to_write Matrix to be written
 * @param root Group to write to
 * @param name Name of the target Attribute
 */
void WriteMatrixX(MatrixX &matrix_to_write, HighFive::Group &root, std::string name);

/**
 * @brief Reads \p matrix_to_read_into from \p root from the attribute with the name \p name
 * @param matrix_to_read_into Matrix to read into
 * @param root Group to read from
 * @param name Name of the target Attribute
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

HighFive::DataSet add_to_expandable_dataset(HighFive::Group &root, const std::string &sub_name,
                                            const std::vector<size_t> &dims, std::vector<size_t> &offset,
                                            bool override = false);

#endif //BACHELOR_THESIS_MYTYPES_H
