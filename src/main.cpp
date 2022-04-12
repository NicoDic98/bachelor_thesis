/**
  * @file       main.cpp
  * @brief      Starting point of the program
  * @author     nico
  * @version    0.0.4
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
#include <iomanip>


/**
 * @brief Tests the Leap Frog integrator
 */
[[maybe_unused]] void test_leap_frog() {
    const int grid_size = 5;
    const int dim = 2;
    const int lambda = int_pow(grid_size, dim);
    const double C{0.1};
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
    const double C{0.1};
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
    for (double inverse_beta = 0.3; inverse_beta < 4.05; inverse_beta += 0.1) {
        test.set_beta(1. / inverse_beta);
        std::cout << "Acceptance rate:" << HMCTest.generate_ensembles(phi0, 20000, 1000, false) << std::endl;
        auto mag = HMCTest.compute_observable(&BaseModel<VectorX>::get_magnetization);
        double m{0.};
        for (auto elem: mag) {
            m += elem;
        }
        m /= static_cast<double>(mag.size());
        std::cout << "Inverse Beta: " << inverse_beta << "\t Magnetization:" << m << std::endl;
        output << inverse_beta << '\t' << m << '\n';
    }
    output.close();

}


void test_multi_level_hmc() {
    const int grid_size = 8;
    const int dim = 2;
    const int lambda = int_pow(grid_size, dim);
    const double C{0.1};
    const double beta{1. / 0.8};

    VectorX phi0(lambda);
    phi0.setZero();
    VectorX h0(lambda);
    h0.setZero();
    VectorX eta0(lambda);
    eta0.setZero();
    std::default_random_engine myengine{42L};

    IsingModel test(beta, h0, eta0, C, dim, 1, grid_size);
    //MultiLevelHMCGenerator mygen(test, {1, 2, 3}, {1, 2, 3}, 2, InterpolationType::Checkerboard, {8, 8, 8},
    //                           {1. / 8, 1. / 8, 1. / 8},
    //                       myengine);



    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d_%m_%Y__%H_%M_%S_");
    std::string my_time{oss.str()};

    for (double inverse_beta = 0.3; inverse_beta < 4.05; inverse_beta += 0.1) {
        test.set_beta(1. / inverse_beta);
        MultiLevelHMCGenerator mygen(test, {1, 2, 3}, {1, 2, 3}, 2, InterpolationType::Checkerboard, {8, 12, 16},
                                     {1. / 8, 1. / 12, 1. / 16},
                                     myengine);
        std::vector<double> acceptance_rates = mygen.generate_ensembles(phi0, 10000, 1000);
        for (auto acceptance_rate: acceptance_rates) {
            std::cout << "Acceptance rate:" << acceptance_rate << std::endl;
        }
        std::string filename{std::string(DATA_DIR).append(my_time).append(std::to_string(inverse_beta)).append(".h5")};
        HighFive::File file(filename,
                            HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        mygen.dumpToH5(file);
        //mygen.dump_observable(&BaseModel<VectorX>::get_magnetization, "magnetization", file);

    }


}

void MultiLevelTime() {
    const int grid_size = 8;
    const int dim = 2;
    const int lambda = int_pow(grid_size, dim);
    const double C{0.1};
    const double beta{1. / 0.8};

    VectorX phi0(lambda);
    phi0.setZero();
    VectorX h0(lambda);
    h0.setZero();
    VectorX eta0(lambda);
    eta0.setZero();
    std::default_random_engine myengine{42L};

    IsingModel test(beta, h0, eta0, C, dim, 1, grid_size);


    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d_%m_%Y__%H_%M_%S_");
    std::string my_time{oss.str()};
    auto my_test = test.get_coarser_model(InterpolationType::Checkerboard);
    my_test->print_interpolation_matrix();
    my_test->print_connectivity_matrix();

    MultiLevelHMCGenerator mygen(test, {1, 2, 3}, {1, 2, 3}, 2, InterpolationType::Checkerboard, {8, 12, 16},
                                 {1. / 8, 1. / 12, 1. / 16},
                                 myengine);
    std::vector<double> acceptance_rates = mygen.generate_ensembles(phi0, 10000, 1000);
    for (auto acceptance_rate: acceptance_rates) {
        std::cout << "Acceptance rate:" << acceptance_rate << std::endl;
    }
    std::string filename{std::string(DATA_DIR).append(my_time).append(std::to_string(1. / beta)).append(".h5")};
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    mygen.dumpToH5(file);
}

void test_hmc_measurements() {
    double inverse_beta{0.8};
    std::string my_time{"12_04_2022__00_34_44_"};
    std::string filename{std::string(DATA_DIR).append(my_time).append(std::to_string(inverse_beta)).append(".h5")};
    HighFive::File file(filename, HighFive::File::ReadOnly);
    auto helper = file.getGroup("level0");//todo see if this step can be removed to be needed
    IsingModel test(helper);

    std::default_random_engine myengine{42L};
    MultiLevelHMCGenerator mygen(test, file, myengine);
    HighFive::File out_file(std::string(DATA_DIR).append("out.h5"),
                            HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    mygen.dumpToH5(out_file);
    mygen.dump_observable(&BaseModel<VectorX>::get_magnetization, "magnetization", out_file);
    mygen.dump_observable(&BaseModel<VectorX>::get_magnetization_squared, "magnetization_squared", out_file);
}

/**
 * @brief Main function
 * @param argc
 * @param argv
 * @return Exit status
 */
int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
    //test_leap_frog();
    //test_HMC(std::string(DATA_DIR).append("HMCTest1.dat"));
    //test_multi_level_hmc();
    test_hmc_measurements();
    //MultiLevelTime();
}


