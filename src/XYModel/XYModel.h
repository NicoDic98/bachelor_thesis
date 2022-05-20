/**
 * @file       XYModel.h
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       18.05.22
 */


#ifndef MULTILEVELHMC_XYMODEL_H
#define MULTILEVELHMC_XYMODEL_H

#include <BaseModel.h>
#include <vector>
#include <MyTypes.h>


class XYModel : public BaseModel<MultiVectorX> {
public:
    XYModel(double beta_, MultiVectorX h_, int dimension_, int neighbour_extent_, int grid_size_);
    /**
     * @brief Coarsening constructor of XYModel
     * @param NewModel Finer model
     * @param InterpolationType_ Interpolation type to use for the coarsening
     */
    XYModel(const XYModel &NewModel, InterpolationType InterpolationType_);

    /**
     * @brief Copy constructor of XYModel
     * @param NewModel Model to copy
     */
    XYModel(const XYModel &NewModel);

    /**
     * @brief Calculates the action for the given field phi
     * @param phi Field
     * @return S(phi) (action)
     */
    double get_action(const MultiVectorX &phi) override;

protected:

    /**
     * @brief Calculates the artificial force for the given field phi
     * @param phi Field
     * @return Artificial force according to artificial hamiltonian
     */
    MultiVectorX get_force(const MultiVectorX &phi) override;

private:
    /**
     * @brief Dimension of the hypercube
     */
    int dimension;
    /**
     * @brief String to be used as key for \a dimension in H5 files
     */
    static const char *dimension_name;

    /**
     * @brief Neighbour extent
     */
    int neighbour_extent;
    /**
     * @brief String to be used as key for \a neighbour_extent in H5 files
     */
    static const char *neighbour_extent_name;

    /**
     * @brief External field
     */
    VectorX h;
    /**
     * @brief String to be used as key for \a h in H5 files
     */
    static const char *h_name;

    /**
     * @brief Symmetric connectivity matrix
     */
    MatrixX k_sym;
    /**
     * @brief String to be used as key for \a k_sym in H5 files
     */
    static const char *k_sym_name;
};


#endif //MULTILEVELHMC_XYMODEL_H
