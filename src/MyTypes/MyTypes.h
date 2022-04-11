/**
 * @file       MyTypes.h
 * @brief      Definition of used types
 * @author     nico
 * @version    0.0.1
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
typedef Eigen::MatrixXd MatrixX;

/**
 * @brief Interpolation types
 */
enum class InterpolationType {
    Checkerboard, Black_White
};


void WriteVectorX(VectorX &vector_to_write, HighFive::Group &root, std::string name);

void ReadVectorX(VectorX &vector_to_read_into, HighFive::Group &root, std::string name);

void WriteMatrixX(MatrixX &matrix_to_write, HighFive::Group &root, std::string name);

void ReadMatrixX(MatrixX &matrix_to_read_into, HighFive::Group &root, std::string name);

#endif //BACHELOR_THESIS_MYTYPES_H
