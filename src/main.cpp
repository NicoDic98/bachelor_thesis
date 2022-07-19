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
#include <XYModel.h>
#include <LeapFrogIntegrator.h>
#include <HMCGenerator.h>
#include <MultiLevelHMCGenerator.h>
#include <Analyzer.h>
#include <fstream>
#include <iomanip>
#include <filesystem>

#include <omp.h>

/**
 * @brief Tests the Leap Frog integrator for the Ising model
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
    //test.print_connectivity_matrix();
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
 * @param filename File to print to
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
        std::cout << "Acceptance configs:" << HMCTest.generate_ensembles(phi0, 20000, 1000, false, 100) << std::endl;
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

/**
 * @brief Tests the basic functionality of MLHMC
 */
void test_multi_level_hmc() {
    const int grid_size = 16;
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
        MultiLevelHMCGenerator mygen(test, {0, 1}, {1, 1}, {-1, -1}, 1, InterpolationType::Checkerboard,
                                     {8, 8},
                                     {1. / 8, 1. / 8},
                                     myengine);
        std::vector<double> acceptance_rates = mygen.generate_ensembles(phi0, 100000, 0);
        for (auto acceptance_rate: acceptance_rates) {
            std::cout << "Acceptance rate:" << acceptance_rate << std::endl;
        }
        std::string filename{std::string(DATA_DIR).append(my_time).append(std::to_string(inverse_beta)).append(".h5")};
        HighFive::File file(filename,
                            HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        mygen.dumpToH5(file);

    }


}

/**
 * @brief Performs measurements and analysis on the given MLHMC generator
 * @tparam configuration_type Gets inferred from \p Gen
 * @param Gen MLHMC generator, which has the ensembles stored
 * @param out_filename Output filename
 * @param remeasure Rather to perform measurements or not
 */
template<class configuration_type>
void DoMultiLevelMeasurements(MultiLevelHMCGenerator<configuration_type> &Gen, const std::string &out_filename,
                              bool remeasure) {
    if (remeasure) {
        HighFive::File out_file(out_filename,
                                HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
        Gen.dump_observable(&BaseModel<configuration_type>::get_magnetization, "magnetization", out_file);
        Gen.dump_observable(&BaseModel<configuration_type>::get_field_squared, "field_squared", out_file);
        Gen.dump_observable(&BaseModel<configuration_type>::get_magnetization_squared, "magnetization_squared",
                            out_file);
        Gen.dump_observable(&BaseModel<configuration_type>::get_energy, "energy", out_file);
        Gen.dump_observable(&BaseModel<configuration_type>::get_energy_squared, "energy_squared", out_file);
        Gen.dump_observable(&BaseModel<configuration_type>::get_vector_length_squared, "vector_length_squared",
                            out_file);
    }
    HighFive::File out_file(out_filename, HighFive::File::ReadWrite);
    int my_size = 1e5;
    int amount_of_sub_samples = 20;
    for (int k = my_size / amount_of_sub_samples; k < my_size; k += my_size / amount_of_sub_samples) {
        Gen.analyze_dataset("magnetization", out_file, -1, -1, k, -1, 200, 20000);
        Gen.analyze_dataset("magnetization_squared", out_file, -1, -1, k, -1, 200, 20000);
        Gen.analyze_dataset("energy", out_file, -1, -1, k, -1, 200, 20000);
        Gen.analyze_dataset("energy_squared", out_file, -1, -1, k, -1, 200, 20000);
    }
    Gen.analyze_dataset("magnetization", out_file, -1, -1, -1, -1, 200, 20000);
    Gen.analyze_dataset("magnetization_squared", out_file, -1, -1, -1, -1, 200, 20000);
    Gen.analyze_dataset("energy", out_file, -1, -1, -1, -1, 200, 20000);
    Gen.analyze_dataset("energy_squared", out_file, -1, -1, -1, -1, 200, 20000);
    //Gen.analyze_dataset("vector_length_squared", out_file, -1, -1, 3000, 200, 400);
}

/**
 * @brief Performs measurements by reading in a given file with name \p filename
 * @param filename Name of file to read in
 * @param model Name of model under which the file was generated
 * @param remeasure Rather to perform measurements or not
 */
void DoMultiLevelMeasurementsFromFile(std::string filename, const std::string &model, bool remeasure) {

    if (!filename.starts_with('/')) {
        filename = std::string(DATA_DIR).append(filename);
    }
    HighFive::File file(filename, HighFive::File::ReadOnly);
    auto helper = file.getGroup("level0");//todo see if this step can be removed to be needed
    if (model == "xy") {
        XYModel test(helper);

        std::default_random_engine myengine{42L};
        MultiLevelHMCGenerator mygen(test, file, myengine);

        filename.insert(filename.rfind('/') + 1, "out_");
        std::string out_filename{filename};
        DoMultiLevelMeasurements(mygen, out_filename, remeasure);
    } else if (model == "ising") {
        IsingModel test(helper);

        std::default_random_engine myengine{42L};
        MultiLevelHMCGenerator mygen(test, file, myengine);

        filename.insert(filename.rfind('/') + 1, "out_");
        std::string out_filename{filename};
        DoMultiLevelMeasurements(mygen, out_filename, remeasure);
    }
}

/**
 * @brief Performs measurements by reading in all h5 files within the directory \p dirname
 * @param dirname Name of the directory to read in
 * @param model Name of model under which the files were generated
 * @param remeasure Rather to perform measurements or not
 */
void DoMultiLevelMeasurementsFromDir(const std::string &dirname, const std::string &model, bool remeasure) {
    for (const auto &file: std::filesystem::directory_iterator(std::string(DATA_DIR).append(dirname))) {
        if (std::string(file.path().filename()).starts_with("out")) {
            continue;
        }
        if (std::string(file.path().filename()).ends_with(".h5")) {
            std::cout << file.path() << std::endl;
            DoMultiLevelMeasurementsFromFile(file.path(), model, remeasure);
        }
    }
}

/**
 * @brief Perform simulation of the Isingmodel at the critical point
 * @param grid_size Grid size on which to simulate (always 2 dimensional)
 * @param nu_pre Amount of pre steps to be performed on the different levels of MLHMC
 * @param nu_post Amount of post steps to be performed on the different levels of MLHMC
 * @param erg_jump_dists Interval of ergodicity jumps on each level of MLHMC, -1 means no jumps
 * @param gamma Amount of gamma repetitions
 * @param int_type Interpolation type
 * @param amount_of_steps Amount of molecular dynamics steps at each level of MLHMC
 * @param step_sizes Step size used in the LeapFrogIntegrator at each level of MLHMC
 * @param id Id to be used for the saved file
 * @param config_amount Amount of configurations to be produced
 * @param dim Dimension of the lattice
 */
void MultiLevelCriticalSimulation(const int grid_size = 16,
                                  std::vector<size_t> nu_pre = {0},
                                  std::vector<size_t> nu_post = {1},
                                  std::vector<int> erg_jump_dists = {-1},
                                  size_t gamma = 1,
                                  InterpolationType int_type = InterpolationType::Checkerboard,
                                  const std::vector<size_t> &amount_of_steps = {6},
                                  const std::vector<double> &step_sizes = {1. / 6.},
                                  size_t id = 0, int config_amount = 30000, int dim = 2) {
    const int lambda = int_pow(grid_size, dim);
    const double C{0.1};
    double beta{4.2};
    if (dim == 2) {
        beta = 0.440686793509772;
    } else if (dim == 3) {
        beta = 0.22165468;
    } else {
        std::cerr << "I do not know the critical point in " << dim << " dimensions.\n";
        exit(dim);
    }
    std::string filename{std::string(DATA_DIR)};

    VectorX phi0(lambda);
    phi0.setZero();
    VectorX h0(lambda);
    h0.setZero();
    VectorX eta0(lambda);
    eta0.setZero();
    std::default_random_engine myengine{42L};

    IsingModel test(beta, h0, eta0, C, dim, 1, grid_size);


    MultiLevelHMCGenerator mygen(test, nu_pre, nu_post, erg_jump_dists, gamma, int_type,
                                 amount_of_steps, step_sizes, myengine);
    std::vector<double> acceptance_rates = mygen.generate_ensembles(phi0, config_amount, 0);
    for (auto acceptance_rate: acceptance_rates) {
        std::cout << "Acceptance rate:" << acceptance_rate << std::endl;
    }

    if (int_type == InterpolationType::Black_White) {
        filename.append("gs").append(std::to_string(grid_size)).append("_")
                .append("Black_White").append("_")
                .append("ga").append(std::to_string(gamma)).append("_")
                .append("l").append(std::to_string(nu_pre.size())).append("_")
                .append("id").append(std::to_string(id)).append(".h5");
    } else if (int_type == InterpolationType::Checkerboard) {
        filename.append("gs").append(std::to_string(grid_size)).append("_")
                .append("Checkerboard").append("_")
                .append("ga").append(std::to_string(gamma)).append("_")
                .append("l").append(std::to_string(nu_pre.size())).append("_")
                .append("id").append(std::to_string(id)).append(".h5");
    }
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    mygen.dumpToH5(file);
}

/**
 * @brief Simulation of the Ising model at the critical point using HMC as a special case of MLHMC
 * @param grid_size Grid size on which to simulate (always 2 dimensional)
 * @param amount_of_steps Amount of molecular dynamics steps at each level of MLHMC
 * @param step_sizes Step size used in the LeapFrogIntegrator at each level of MLHMC
 * @param id Id to be used for the saved file
 * @param config_amount Amount of configurations to be produced
 * @param dim Dimension of the lattice
 */
void HMCCriticalSimulation(int grid_size = 16, const size_t &amount_of_steps = 6,
                           const double step_sizes = 1. / 6., size_t id = 0, int config_amount = 30000, int dim = 2) {
    MultiLevelCriticalSimulation(grid_size, {0}, {1},
                                 {-1}, 1, InterpolationType::Checkerboard,
                                 {amount_of_steps},
                                 {step_sizes}, id, config_amount, dim);
}

/**
 * @brief Tests the Leap Frog integrator for the XY model
 */
void test_leap_frog_XY() {
    const int grid_size = 5;

    const int dim = 2;
    const int lambda = int_pow(grid_size, dim);
    const double beta{1.27};
    std::string filename{std::string(DATA_DIR)};

    VectorX temp(lambda);
    temp.setOnes();
    MultiVectorX phi0;
    //phi0.push_back(cos(temp.array()));
    //phi0.push_back(sin(temp.array()));
    phi0.push_back(temp);
    temp.setZero();
    phi0.push_back(temp);

    temp.setZero();
    MultiVectorX h0;
    h0.push_back(temp);
    h0.push_back(temp);
    temp.setRandom();
    MultiVectorX pi0;
    pi0.push_back(temp);
    pi0.push_back(temp);
    std::default_random_engine myengine{42L};

    XYModel test(beta, -6, h0, dim, 1, grid_size);

    LeapFrogIntegrator leapTest(test);

    double S_start = test.get_artificial_energy(phi0, pi0);
    std::vector<double> S_error;
    std::vector<int> MD;
    for (int i = 10; i < 200; i += 20) {
        MD.push_back(i);
        MultiVectorX phiNew;
        MultiVectorX piNew;
        for (const auto &elem: phi0) {
            phiNew.push_back(elem);
        }
        for (const auto &elem: pi0) {
            piNew.push_back(elem);
        }


        leapTest.integrate(i, 1. / i, phiNew, piNew);
        double S_end = test.get_artificial_energy(phiNew, piNew);
        S_error.push_back(abs((S_start - S_end) / S_start));
        //std::cout << S_start<<std::endl;
        //std::cout << S_end<<std::endl;
        std::cout << i << "\t: " << abs((S_start - S_end) / S_start) << std::endl;
    }
}

/**
 * @brief Perform simulation of the XY model at inverse temperature \p beta
 * @param grid_size Grid size on which to simulate (always 2 dimensional)
 * @param beta Inverse temperature at which to simulate
 * @param nu_pre Amount of pre steps to be performed on the different levels of MLHMC
 * @param nu_post Amount of post steps to be performed on the different levels of MLHMC
 * @param erg_jump_dists Interval of ergodicity jumps on each level of MLHMC, -1 means no jumps
 * @param gamma Amount of gamma repetitions
 * @param eta Self introduced lagrange multiplier, which didn't yield promising results
 * @param int_type Interpolation type
 * @param amount_of_steps Amount of molecular dynamics steps at each level of MLHMC
 * @param step_sizes Step size used in the LeapFrogIntegrator at each level of MLHMC
 * @param id Id to be used for the saved file
 */
void MultiLevelSimulationXY(const int grid_size = 16,
                            const double beta = 2,
                            std::vector<size_t> nu_pre = {0},
                            std::vector<size_t> nu_post = {1},
                            std::vector<int> erg_jump_dists = {-1},
                            size_t gamma = 1,
                            double eta = -2,
                            InterpolationType int_type = InterpolationType::Checkerboard,
                            const std::vector<size_t> &amount_of_steps = {6},
                            const std::vector<double> &step_sizes = {1. / 6.},
                            size_t id = 0) {
    const int dim = 2;
    const int lambda = int_pow(grid_size, dim);
    std::string filename{std::string(DATA_DIR)};

    VectorX temp(lambda);
    temp.setOnes();
    MultiVectorX phi0;
    //phi0.push_back(cos(temp.array()));
    //phi0.push_back(sin(temp.array()));
    phi0.push_back(temp);
    temp.setZero();
    phi0.push_back(temp);

    temp.setZero();
    MultiVectorX h0;
    h0.push_back(temp);
    h0.push_back(temp);
    std::default_random_engine myengine{42L};

    XYModel test(beta, eta, h0, dim, 1, grid_size);


    MultiLevelHMCGenerator mygen(test, nu_pre, nu_post, erg_jump_dists, gamma, int_type,
                                 amount_of_steps, step_sizes, myengine);
    std::vector<double> acceptance_rates = mygen.generate_ensembles(phi0, 10000, 10000);
    for (auto acceptance_rate: acceptance_rates) {
        std::cout << "Acceptance rate:" << acceptance_rate << std::endl;
    }

    if (int_type == InterpolationType::Black_White) {
        filename.append("gs").append(std::to_string(grid_size)).append("_")
                .append("Black_White").append("_")
                .append("ga").append(std::to_string(gamma)).append("_")
                .append("l").append(std::to_string(nu_pre.size())).append("_")
                .append("id").append(std::to_string(id)).append(".h5");
    } else if (int_type == InterpolationType::Checkerboard) {
        filename.append("gs").append(std::to_string(grid_size)).append("_")
                .append("Checkerboard").append("_")
                .append("ga").append(std::to_string(gamma)).append("_")
                .append("l").append(std::to_string(nu_pre.size())).append("_")
                .append("id").append(std::to_string(id)).append(".h5");
    }
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    mygen.dumpToH5(file);
}

/**
 * @brief Performs simulation of the XY model over a range of temperatures
 * @param grid_size Grid size on which to simulate (always 2 dimensional)
 * @param amount_of_steps Amount of molecular dynamics steps at each level of MLHMC
 * @param step_sizes Step size used in the LeapFrogIntegrator at each level of MLHMC
 */
void HMCCriticalSimulationXY(int grid_size = 16, const size_t &amount_of_steps = 6,
                             const double step_sizes = 1. / 6.) {

    double inverse_beta = 0.3;
    size_t id = 0;
    while (inverse_beta < 2) {
        MultiLevelSimulationXY(grid_size, 1. / inverse_beta, {0}, {1},
                               {-1}, 1, -3. / (sqrt(inverse_beta)), InterpolationType::Checkerboard,
                               {amount_of_steps},
                               {step_sizes}, id);
        id++;
        inverse_beta += 0.1;
    }


}

/**
 * @brief Main function
 * @param argc Argument count
 * @param argv Arguments
 * @return Exit status
 */
int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
    if (argc >= 2) {
        for (int l = 1; l < argc; ++l) {
            if (std::string(argv[l]) == "license") {
                std::string line;
                std::ifstream license_file(std::string(PROJ_DIR).append("LICENSE"));
                if (license_file.is_open()) {
                    while (license_file) {
                        std::getline(license_file, line);
                        std::cout << line << '\n';
                    }
                }
            }
        }
    }
    //test_leap_frog();
    //test_HMC(std::string(DATA_DIR).append("HMCTest1.dat"));
    //test_multi_level_hmc();
    //test_hmc_measurements();
    //DoMultiLevelMeasurementsFromDir(std::string("3d_new_c_0_3"), std::string("ising"), true);
    size_t i{42};
    //HMCCriticalSimulation(2, 3, 1. / 3., i++, 100000, 3);
    //HMCCriticalSimulation(4, 4, 1. / 4., i++, 100000, 3);
    //HMCCriticalSimulation(8, 8, 1. / 8., i++, 100000, 3);
    //HMCCriticalSimulationXY(16, 12, 1. / 12.);
    std::vector<size_t> nu_pre = {0, 0, 0, 0};
    std::vector<size_t> nu_post = {1, 1, 1, 1};
    std::vector<int> erg_jump_dists = {-1, -1, -1, -1};
    std::vector<size_t> amount_of_steps = {10, 10, 10, 10};
    std::vector<double> step_sizes = {1. / 10., 1. / 10., 1. / 10., 1. / 10.};
    //for (size_t l = 1; l < 17; l *= 4) {
    //nu_pre.push_back(1);
    //nu_post.push_back(1);
    //erg_jump_dists.push_back(-1);
    //amount_of_steps.push_back(l * 3);
    //step_sizes.push_back(1. / (static_cast<double>(l) * 3.));
    MultiLevelCriticalSimulation(32, nu_pre, nu_post,
                                 erg_jump_dists, 2, InterpolationType::Checkerboard,
                                 amount_of_steps,
                                 step_sizes, i++, 100000, 2);
    nu_pre = {0, 0, 0};
    nu_post = {1, 1, 1};
    erg_jump_dists = {-1, -1, -1};
    amount_of_steps = {8, 12, 18};
    step_sizes = {1. / 8., 1. / 12., 1. / 18.};
    MultiLevelCriticalSimulation(8, nu_pre, nu_post,
                                 erg_jump_dists, 1, InterpolationType::Checkerboard,
                                 amount_of_steps,
                                 step_sizes, i++, 100000, 3);
    //}*/
}


