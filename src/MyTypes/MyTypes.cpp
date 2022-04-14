/**
 * @file       MyTypes.cpp
 * @brief
 * @author     nico
 * @version    0.0.1
 * @date       11.04.22
 */
#include "MyTypes.h"

void WriteVectorX(VectorX &vector_to_write, HighFive::Group &root, std::string name) {
    std::vector<size_t> offset;
    assert(vector_to_write.rows() >= 0);

    HighFive::DataSet target_dataset = add_to_expandable_dataset(
            root, name,
            {static_cast<unsigned long>(vector_to_write.rows())},
            offset, true);

    std::vector<size_t> count{static_cast<unsigned long>(vector_to_write.rows())};
    target_dataset.select(offset, count).write_raw(vector_to_write.data());
}

void ReadVectorX(VectorX &vector_to_read_into, HighFive::Group &root, std::string name) {
    auto temp = root.getDataSet(name);
    std::vector<double> my_buffer;
    temp.read(my_buffer);
    vector_to_read_into.resize(my_buffer.size());
    vector_to_read_into = VectorX::Map(&my_buffer[0], my_buffer.size());
}

void WriteMatrixX(MatrixX &matrix_to_write, HighFive::Group &root, std::string name) {
    std::vector<size_t> offset;
    assert(matrix_to_write.rows() >= 0);
    assert(matrix_to_write.cols() >= 0);

    HighFive::DataSet target_dataset = add_to_expandable_dataset(
            root, name,
            {static_cast<unsigned long>(matrix_to_write.rows()), static_cast<unsigned long>(matrix_to_write.cols())},
            offset, true, {0, 1});

    std::vector<size_t> count{static_cast<unsigned long>(matrix_to_write.rows()),
                              static_cast<unsigned long>(matrix_to_write.cols())};
    target_dataset.select(offset, count).write_raw(matrix_to_write.data());
}

void ReadMatrixX(MatrixX &matrix_to_read_into, HighFive::Group &root, std::string name) {
    auto temp = root.getDataSet(name);
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
                                            bool override, const std::vector<size_t> &dims_to_resize) {
    offset.clear();
    offset.resize(dims.size());
    std::fill(offset.begin(), offset.end(), 0);

    if (root.exist(sub_name)) {
        HighFive::DataSet target_dataset = root.getDataSet(sub_name);
        std::vector<size_t> shape{target_dataset.getDimensions()};

        if (shape.size() == dims.size()) {
            for (int i = 0; i < shape.size(); ++i) {
                for (auto resize_dim: dims_to_resize) {
                    if (i == resize_dim) {
                        continue;
                    }
                }
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
        for (auto resize_dim: dims_to_resize) {
            max_dims[resize_dim] = HighFive::DataSpace::UNLIMITED;
        }
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
