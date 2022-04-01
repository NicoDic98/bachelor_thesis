/**
  * @file       main.cpp
  * @brief      Starting point of the program
  * @author     nico
  * @version    0.0.1
  * @date       22.03.22
  */
#include <iostream>
#include <random>
#include <SrcConfig.h>
#include <IsingModel.h>
#include <LeapFrogIntegrator.h>
#include <HMCGenerator.h>
#include <MultiLevelHMCGenerator.h>
#include <fstream>
#include <highfive/H5Easy.hpp>

/**
 * @brief Tests the Leap Frog integrator
 */
void test_leap_frog() {
    const int grid_size = 5;
    const int dim = 2;
    const int lambda = int_pow(grid_size, dim);
    const double C{4.1};
    const double beta{1.0};

    VectorX phi0(lambda);
    phi0.setRandom();
    VectorX pi0(lambda);
    pi0.setRandom();
    VectorX h0(lambda);
    h0.setRandom();
    VectorX eta0(lambda);
    eta0.setRandom();

    IsingModel test(beta, h0, eta0, C, dim, 1, grid_size);

    LeapFrogIntegrator leapTest(test);
    test.print_connectivity_matrix();
    int numMD = 10;
    phi0(0) = 1;
    phi0(1) = 2;
    phi0(2) = -0.1;
    phi0(3) = -1.5;

    double S_start = test.get_action(phi0) + pi0.dot(pi0) / 2.;
    std::vector<double> S_error;
    std::vector<int> MD;
    for (int i = 10; i < 200; i += 20) {
        MD.push_back(i);
        VectorX phiNew(phi0);
        VectorX piNew(pi0);

        leapTest.integrate(i, 1. / i, phiNew, piNew);
        double S_end = test.get_action(phiNew) + piNew.dot(piNew) / 2.;
        S_error.push_back(abs((S_start - S_end) / S_start));
        //std::cout << S_start<<std::endl;
        //std::cout << S_end<<std::endl;
        std::cout << i << "\t: " << abs((S_start - S_end) / S_start) << std::endl;
    }
}

/**
 * @brief Tests the HMC, and prints some results to filename
 * @param filename
 */
void test_HMC(const std::string &filename) {
    const int grid_size = 5;
    const int dim = 2;
    const int lambda = int_pow(grid_size, dim);
    const double C{4.1};
    const double beta{0.5};

    VectorX phi0(lambda);
    phi0.setRandom();
    VectorX h0(lambda);
    h0.setZero();
    VectorX eta0(lambda);
    eta0.setZero();

    IsingModel test(beta, h0, eta0, C, dim, 1, grid_size);

    std::default_random_engine myengine{42L};
    HMCGenerator HMCTest(test, 8, 1. / 8, myengine);

    std::ofstream output(filename);
    if (!output) {
        std::cerr << filename << " can't be opened!\n";
        exit(-42);
    }
    for (double inverse_beta = 0.5; inverse_beta < 1.05; inverse_beta += 0.1) {
        test.set_beta(1. / inverse_beta);
        std::cout << "Acceptance rate:" << HMCTest.generate_ensembles(phi0, 20000, 1000) << std::endl;
        double m = HMCTest.compute_magnetization();
        std::cout << "Beta: " << 1. / inverse_beta << "\t Magnetization:" << m << std::endl;
        output << inverse_beta << '\t' << m << '\n';
    }
    output.close();

}


void test_multi_level_hmc() {
    const int grid_size = 8;
    const int dim = 2;
    const int lambda = int_pow(grid_size, dim);
    const double C{4.1};
    const double beta{0.5};

    VectorX h0(lambda);
    h0.setZero();
    VectorX eta0(lambda);
    eta0.setZero();
    std::default_random_engine myengine{42L};

    IsingModel test(beta, h0, eta0, C, dim, 1, grid_size);
    MultiLevelHMCGenerator(test, {1, 2, 3}, {1, 2, 3}, 1, InterpolationType::Checkerboard, {8, 8, 8},
                           {1. / 8, 1. / 8, 1. / 8},
                           myengine);
}

/**
 * @brief Main function
 * @param argc
 * @param argv
 * @return Exit status
 */
int main(int argc, char *argv[]) {
    /* std::cout << "Hello\n";
     using namespace HighFive;
 // we create a new hdf5 file
     File file(std::string(DATA_DIR).append("new_file.h5"), File::ReadWrite | File::Create | File::Truncate);

     std::vector<int> data(50, 1);

 // let's create a dataset of native integer with the size of the vector 'data'
     DataSet dataset = file.createDataSet<int>("/dataset_one",  DataSpace::From(data));

 // let's write our vector of int to the HDF5 dataset
     dataset.write(data);

 // read back
     std::vector<int> result;
     dataset.read(result);*/
    //test_leap_frog();
    //test_HMC(std::string(DATA_DIR).append("Test.dat"));
    test_multi_level_hmc();
}
