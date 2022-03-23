/**
  * @file       main.cpp
  * @brief      Starting point of the program
  * @author     nico
  * @version    0.0.1
  * @date       22.03.22
  */
#include <iostream>
#include <SrcConfig.h>
#include <IsingModel.h>


int main(int argc, char *argv[]) {
    std::cout << "Hello\n";
    const int grid_size = 10;
    const int dim = 3;
    const int size = int_pow(grid_size, dim);
    IsingModel test(1., VectorX(size).setOnes(), VectorX(size).setZero(), 0.,
                    dim, 1, grid_size);
    test.print_connectivity_matrix();
    test.get_force(VectorX(1000).setOnes());
}
