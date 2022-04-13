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

HighFive::DataSet add_to_expandable_dataset(HighFive::Group &root, const std::string &sub_name,
                                            const std::vector<size_t> &dims, std::vector<size_t> &offset,
                                            bool override) {
    offset.clear();
    offset.resize(dims.size());
    std::fill(offset.begin(), offset.end(), 0);

    if (root.exist(sub_name)) {
        HighFive::DataSet target_dataset = root.getDataSet(sub_name);
        std::vector<size_t> shape{target_dataset.getDimensions()};

        if (shape.size() == dims.size()) {
            for (int i = 0; i < shape.size(); ++i) {
                if (shape[i] != dims[i]) {
                    std::cerr << "Vector sizes not matching\n";
                    exit(-1);
                }
            }
            if (override) {
                shape[0] = dims[0];
                target_dataset.resize(shape);
            } else {
                offset[0] = shape[0];//extend always in first dimension
                shape[0] += dims[0];
                target_dataset.resize(shape);
            }

        } else {
            std::cerr << "Can't extend dataset\n";
            exit(-1);
        }

        return target_dataset;
    } else {
        auto max_dims(dims);
        max_dims[0] = HighFive::DataSpace::UNLIMITED;
        auto my_dataspace = HighFive::DataSpace(dims, max_dims);

        HighFive::DataSetCreateProps props;

        std::vector<hsize_t> chunk_size(dims.size());
        chunk_size[0] = 1;
        for (int i = 1; i < chunk_size.size(); ++i) {
            chunk_size[i] = dims[i];
        }

        props.add(HighFive::Chunking(chunk_size));
        HighFive::DataSet target_dataset = root.createDataSet(sub_name, my_dataspace,
                                                              HighFive::create_datatype<double>(), props);

        return target_dataset;
    }
}
