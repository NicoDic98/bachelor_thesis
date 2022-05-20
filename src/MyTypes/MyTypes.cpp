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

void WriteMultiVectorX(MultiVectorX &vector_to_write, HighFive::Group &root, std::string name) {
    std::vector<size_t> offset;
    assert(vector_to_write.size() >= 0);
    for (const auto &elem: vector_to_write) {
        assert(elem.rows() >= 0);
    }


    HighFive::DataSet target_dataset = add_to_expandable_dataset(
            root, name,
            {static_cast<unsigned long>(vector_to_write.size(), vector_to_write[0].rows())},
            offset, true);

    std::vector<size_t> count{static_cast<unsigned long>(1, vector_to_write[0].rows())};
    for (int l = 0; l < vector_to_write.size(); ++l) {
        offset[0] = l;
        target_dataset.select(offset, count).write_raw(vector_to_write[l].data());
    }

}

void ReadMultiVectorX(MultiVectorX &vector_to_read_into, HighFive::Group &root, std::string name) {
    auto temp = root.getDataSet(name);
    std::vector<std::vector<double>> my_buffer;
    temp.read(my_buffer);

    vector_to_read_into.clear();
    for (auto & item : my_buffer) {
        vector_to_read_into.push_back(VectorX::Map(&item[0], item.size()));
    }
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
                bool is_resize_dim = false;
                for (auto resize_dim: dims_to_resize) {
                    if (i == resize_dim) {
                        is_resize_dim = true;
                        break;
                    }
                }
                if (is_resize_dim) {
                    continue;
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

void fill_connectivity_matrix(int grid_size, int dimension, int neighbour_extent, MatrixX &k_sym) {
    k_sym.setZero();

    long lambda = k_sym.rows();

    for (long m = 0; m < lambda; ++m) {
        long base_offset{1};
        for (int i = 0; i < dimension; ++i) {
            for (long j = 1; j <= neighbour_extent; ++j) {
                k_sym(m, (m + base_offset * neighbour_extent) % (base_offset * grid_size)
                         + (base_offset * grid_size) * (m / (base_offset * grid_size))) += 1;
                // this additional + is needed in cas of eg ((0,1),(2,3)) to get the correct connections for 2
                k_sym(m, (m - base_offset * neighbour_extent + (base_offset * grid_size)) % (base_offset * grid_size)
                         + (base_offset * grid_size) * (m / (base_offset * grid_size))) += 1;
            }
            base_offset *= grid_size;
        }

    }
}

void fill_interpolation_matrix(InterpolationType InterpolationType_, long fine_size, int dimension,
                               MatrixX &InterpolationMatrix) {
    int coarse_grid_side_length{1};
    long coarse_lambda;
    int fine_grid_side_length = static_cast<int>(pow(static_cast<double>(fine_size), 1. / dimension));
    switch (InterpolationType_) {
        case InterpolationType::Checkerboard:
            if (int_pow(fine_grid_side_length, dimension) != fine_size) {
                if (int_pow(fine_grid_side_length + 1, dimension) == fine_size) {
                    fine_grid_side_length += 1;
                } else if (int_pow(fine_grid_side_length - 1, dimension) == fine_size) {
                    fine_grid_side_length -= 1;
                } else {
                    std::cerr << "Can't determine fine grid side length from fine_size = "
                              << fine_size << " and dimension = " << dimension << '\n';
                    exit(-42);
                }
            }


            /*
             * Compute new sizes
             */
            if (fine_grid_side_length % 2 == 0) {
                coarse_grid_side_length = fine_grid_side_length / 2;
            } else {
                coarse_grid_side_length = (fine_grid_side_length + 1) / 2;
                //TODO: think more about this case
            }

            coarse_lambda = int_pow(coarse_grid_side_length, dimension);

            InterpolationMatrix.resize(fine_size, coarse_lambda);
            InterpolationMatrix.setZero();


            for (long m = 0; m < InterpolationMatrix.rows(); ++m) {
                long base_offset{1};
                long coarse_offset;
                double factor{1};
                /*
                 * Get ids of fine level sources of the current m
                 */
                std::vector<long> fine_interpolation_ids{m};
                for (int i = 0; i < dimension; ++i) {
                    if (((m / base_offset) % fine_grid_side_length) % 2 == 1) {
                        factor *= 0.5;
                        std::vector<long> temp{fine_interpolation_ids};
                        fine_interpolation_ids.clear();
                        for (long elem: temp) {
                            fine_interpolation_ids.push_back(
                                    (elem + base_offset) % (base_offset * fine_grid_side_length) +
                                    (base_offset * fine_grid_side_length) *
                                    (elem / (base_offset * fine_grid_side_length)));
                            fine_interpolation_ids.push_back(
                                    (elem - base_offset + (base_offset * fine_grid_side_length)) %
                                    (base_offset * fine_grid_side_length) +
                                    (base_offset * fine_grid_side_length) *
                                    (elem / (base_offset * fine_grid_side_length)));
                        }
                    }

                    base_offset *= fine_grid_side_length;
                }
                //assert(factor == fine_interpolation_ids.size());

                /*
                 * Fill interpolation matrix, by calculating the coarse level id of each source
                 */
                for (auto elem: fine_interpolation_ids) {
                    base_offset = 1;
                    coarse_offset = 1;
                    long coarse_id{0};
                    for (int i = 0; i < dimension; ++i) {
                        coarse_id +=
                                ((elem % (base_offset * fine_grid_side_length)) / (base_offset * 2)) * coarse_offset;

                        base_offset *= fine_grid_side_length;
                        coarse_offset *= coarse_grid_side_length;
                    }
                    InterpolationMatrix(m, coarse_id) = factor;
                }

            }
            break;
        case InterpolationType::Black_White:
            /*
             * Check for d=2
             */
            if (dimension != 2) {
                std::cerr << "Black_White only suppoerted for d=2\n";
                exit(-1);
            }

            /**
             * Compute new sizes
             */
            if (fine_size % 2 != 0) {
                std::cerr << "Black_White only suppoerted for even side lengths\n";
                exit(-1);
            }

            coarse_lambda = fine_size / 2;

            InterpolationMatrix.resize(fine_size, coarse_lambda);
            InterpolationMatrix.setZero();

            bool is_square_lattice = false;

            if (int_pow(fine_grid_side_length, dimension) != fine_size) {
                if (int_pow(fine_grid_side_length + 1, dimension) == fine_size) {
                    fine_grid_side_length += 1;
                    is_square_lattice = true;
                } else if (int_pow(fine_grid_side_length - 1, dimension) == fine_size) {
                    fine_grid_side_length -= 1;
                    is_square_lattice = true;
                }
            } else {
                is_square_lattice = true;
            }

            if (is_square_lattice) {
                //in this case fine_grid_side_length is well-defined
                for (int m = 0; m < InterpolationMatrix.rows(); ++m) {
                    if ((m / fine_grid_side_length) % 2 == 0) {
                        if ((m % fine_grid_side_length) % 2 == 0) {
                            InterpolationMatrix(m,
                                                (((m + 1) % fine_grid_side_length +
                                                  (m / fine_grid_side_length) * fine_grid_side_length) - 1) / 2) = 0.25;
                            InterpolationMatrix(m,
                                                (((m - 1 + fine_grid_side_length) % fine_grid_side_length +
                                                  (m / fine_grid_side_length) * fine_grid_side_length) - 1) / 2) = 0.25;
                            InterpolationMatrix(m, ((m + fine_grid_side_length) % fine_size) / 2) = 0.25;
                            InterpolationMatrix(m, ((m - fine_grid_side_length + fine_size) % fine_size) / 2) = 0.25;
                        } else {
                            InterpolationMatrix(m, (m - 1) / 2) = 1.;
                        }
                    } else {
                        if (m % 2 == 0) {
                            InterpolationMatrix(m, m / 2) = 1.;
                        } else {
                            InterpolationMatrix(m,
                                                ((m + 1) % fine_grid_side_length +
                                                 (m / fine_grid_side_length) * fine_grid_side_length) / 2) = 0.25;
                            InterpolationMatrix(m,
                                                ((m - 1 + fine_grid_side_length) % fine_grid_side_length +
                                                 (m / fine_grid_side_length) * fine_grid_side_length) / 2) = 0.25;
                            InterpolationMatrix(m, (((m + fine_grid_side_length) % fine_size) - 1) / 2) = 0.25;
                            InterpolationMatrix(m,
                                                (((m - fine_grid_side_length + fine_size) % fine_size) - 1) / 2) = 0.25;
                        }
                    }
                }
            } else {
                //in this case fine_grid_side_length is NOT well-defined
                coarse_grid_side_length = static_cast<int>(pow(static_cast<double>(fine_size) / 2, 1. / dimension));
                if (int_pow(coarse_grid_side_length, dimension) != fine_size / 2) {
                    if (int_pow(coarse_grid_side_length + 1, dimension) == fine_size / 2) {
                        coarse_grid_side_length += 1;
                    } else if (int_pow(coarse_grid_side_length - 1, dimension) == fine_size / 2) {
                        coarse_grid_side_length -= 1;
                    } else {
                        std::cerr << "Can't determine coarse grid side length from fine_size = "
                                  << fine_size << " and dimension = " << dimension << '\n';
                        exit(-42);
                    }
                }

                for (int m = 0; m < InterpolationMatrix.rows(); ++m) {
                    if ((m / coarse_grid_side_length) % 2 == 1) {
                        InterpolationMatrix(m, m - (m / coarse_grid_side_length) * coarse_grid_side_length +
                                               (((m / coarse_grid_side_length) - 1) / 2) *
                                               coarse_grid_side_length) = 1.;
                    } else {
                        int neighbourIndex = (m / coarse_grid_side_length - 1) * coarse_grid_side_length +
                                             m % coarse_grid_side_length;
                        neighbourIndex = static_cast<int>((neighbourIndex + fine_size) % fine_size);
                        InterpolationMatrix(m, neighbourIndex -
                                               (neighbourIndex / coarse_grid_side_length) * coarse_grid_side_length +
                                               (((neighbourIndex / coarse_grid_side_length) - 1) / 2) *
                                               coarse_grid_side_length) = 0.25;
                        neighbourIndex = (m / coarse_grid_side_length - 1) * coarse_grid_side_length +
                                         (m + 1) % coarse_grid_side_length;
                        neighbourIndex = static_cast<int>((neighbourIndex + fine_size) % fine_size);
                        InterpolationMatrix(m, neighbourIndex -
                                               (neighbourIndex / coarse_grid_side_length) * coarse_grid_side_length +
                                               (((neighbourIndex / coarse_grid_side_length) - 1) / 2) *
                                               coarse_grid_side_length) = 0.25;
                        neighbourIndex = (m / coarse_grid_side_length + 1) * coarse_grid_side_length +
                                         m % coarse_grid_side_length;
                        InterpolationMatrix(m, neighbourIndex -
                                               (neighbourIndex / coarse_grid_side_length) * coarse_grid_side_length +
                                               (((neighbourIndex / coarse_grid_side_length) - 1) / 2) *
                                               coarse_grid_side_length) = 0.25;
                        neighbourIndex = (m / coarse_grid_side_length + 1) * coarse_grid_side_length +
                                         (m + 1) % coarse_grid_side_length;
                        InterpolationMatrix(m, neighbourIndex -
                                               (neighbourIndex / coarse_grid_side_length) * coarse_grid_side_length +
                                               (((neighbourIndex / coarse_grid_side_length) - 1) / 2) *
                                               coarse_grid_side_length) = 0.25;

                    }
                }
            }
            break;
    }
}
