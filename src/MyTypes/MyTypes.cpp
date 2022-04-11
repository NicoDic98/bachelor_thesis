/**
 * @file       MyTypes.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       11.04.22
 */
#include "MyTypes.h"
void WriteVectorX(VectorX &vector_to_write, HighFive::Group &root, std::string name) {
    if (root.hasAttribute(name)) {
        root.deleteAttribute(name);
    }
    auto temp = root.createAttribute(name, HighFive::DataSpace(vector_to_write.rows()),
                                     HighFive::create_datatype<double>());
    temp.write_raw(vector_to_write.data());
}

void ReadVectorX(VectorX &vector_to_read_into, HighFive::Group &root, std::string name) {
    auto temp = root.getAttribute(name);
    std::vector<double> my_buffer;
    temp.read(my_buffer);

    VectorX vec_temp(my_buffer.size());
    vec_temp = VectorX::Map(&my_buffer[0], my_buffer.size());
    vector_to_read_into = vec_temp;
}

void WriteMatrixX(MatrixX &matrix_to_write, HighFive::Group &root, std::string name) {
    if (root.hasAttribute(name)) {
        root.deleteAttribute(name);
    }
    auto temp = root.createAttribute(name,
                                     HighFive::DataSpace(matrix_to_write.rows(), matrix_to_write.cols()),
                                     HighFive::create_datatype<double>());
    temp.write_raw(matrix_to_write.data());
}

void ReadMatrixX(MatrixX &matrix_to_read_into, HighFive::Group &root, std::string name) {
    auto temp = root.getAttribute(name);
    std::vector<std::vector<double>> my_buffer;
    temp.read(my_buffer);

    MatrixX mat_temp(my_buffer.size(), my_buffer[0].size());
    for (int i = 0; i < my_buffer.size(); ++i) {
        mat_temp.row(i) = VectorX::Map(&my_buffer[i][0], mat_temp.cols());
    }
    matrix_to_read_into = mat_temp;
}
