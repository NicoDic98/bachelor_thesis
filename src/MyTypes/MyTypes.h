/**
 * @file       MyTypes.h
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       24.03.22
 */


#ifndef BACHELOR_THESIS_MYTYPES_H
#define BACHELOR_THESIS_MYTYPES_H

#include <eigen3/Eigen/Dense>

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

#endif //BACHELOR_THESIS_MYTYPES_H
