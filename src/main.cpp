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
#include <Analyzer.h>
#include <fstream>
#include <iomanip>
#include <hip/hip_runtime.h>

#define HIP_ASSERT(x) (assert((x)==hipSuccess))


#define WIDTH     1024
#define HEIGHT    1024

#define NUM       (WIDTH*HEIGHT)

#define THREADS_PER_BLOCK_X  16
#define THREADS_PER_BLOCK_Y  16
#define THREADS_PER_BLOCK_Z  1

__global__ void
vectoradd_float(float *__restrict__ a, const float *__restrict__ b, const float *__restrict__ c, int width,
                int height) {

    int x = hipBlockDim_x * hipBlockIdx_x + hipThreadIdx_x;
    int y = hipBlockDim_y * hipBlockIdx_y + hipThreadIdx_y;

    int i = y * width + x;
    if (i < (width * height)) {
        a[i] = b[i] + c[i];
    }


}

int test_hip() {
    float *hostA;
    float *hostB;
    float *hostC;

    float *deviceA;
    float *deviceB;
    float *deviceC;

    hipDeviceProp_t devProp;
    auto ret = (hipGetDeviceProperties(&devProp, 0));
    HIP_ASSERT(ret);
    std::cout << " System minor " << devProp.minor << std::endl;
    std::cout << " System major " << devProp.major << std::endl;
    std::cout << " agent prop name " << devProp.name << std::endl;


    std::cout << "hip Device prop succeeded " << std::endl;


    int i;
    int errors;

    hostA = (float *) malloc(NUM * sizeof(float));
    hostB = (float *) malloc(NUM * sizeof(float));
    hostC = (float *) malloc(NUM * sizeof(float));

    // initialize the input data
    for (i = 0; i < NUM; i++) {
        hostB[i] = (float) i;
        hostC[i] = (float) i * 100.0f;
    }

    ret = (hipMalloc((void **) &deviceA, NUM * sizeof(float)));
    HIP_ASSERT(ret);
    ret = (hipMalloc((void **) &deviceB, NUM * sizeof(float)));
    HIP_ASSERT(ret);
    ret = (hipMalloc((void **) &deviceC, NUM * sizeof(float)));
    HIP_ASSERT(ret);

    ret = (hipMemcpy(deviceB, hostB, NUM * sizeof(float), hipMemcpyHostToDevice));
    HIP_ASSERT(ret);
    ret = (hipMemcpy(deviceC, hostC, NUM * sizeof(float), hipMemcpyHostToDevice));
    HIP_ASSERT(ret);


    hipLaunchKernelGGL(vectoradd_float,
                       dim3(WIDTH / THREADS_PER_BLOCK_X, HEIGHT / THREADS_PER_BLOCK_Y),
                       dim3(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y),
                       0, 0,
                       deviceA, deviceB, deviceC, WIDTH, HEIGHT);


    ret = (hipMemcpy(hostA, deviceA, NUM * sizeof(float), hipMemcpyDeviceToHost));
    HIP_ASSERT(ret);

    // verify the results
    errors = 0;
    for (i = 0; i < NUM; i++) {
        if (hostA[i] != (hostB[i] + hostC[i])) {
            errors++;
        }
    }
    if (errors != 0) {
        std::cerr << "FAILED: " << errors << " errors";
    } else {
        std::cout << "PASSED!" << std::endl;
    }

    ret = (hipFree(deviceA));
    HIP_ASSERT(ret);
    ret = (hipFree(deviceB));
    HIP_ASSERT(ret);
    ret = (hipFree(deviceC));
    HIP_ASSERT(ret);

    free(hostA);
    free(hostB);
    free(hostC);

    //hipResetDefaultAccelerator();

    return errors;
}

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

    }


}

template<class configuration_type>
void DoMultiLevelMeasurements(MultiLevelHMCGenerator<configuration_type> &Gen, const std::string &out_filename) {
    HighFive::File out_file(out_filename,
                            HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    Gen.dump_observable(&BaseModel<VectorX>::get_magnetization, "magnetization", out_file);
    Gen.analyze_dataset("magnetization", out_file, 100, 200, 200);
    Gen.dump_observable(&BaseModel<VectorX>::get_magnetization_squared, "magnetization_squared", out_file);
    Gen.analyze_dataset("magnetization_squared", out_file, 100, 200, 200);
    Gen.dump_observable(&BaseModel<VectorX>::get_energy, "energy", out_file);
    Gen.analyze_dataset("energy", out_file, 100, 200, 200);
    Gen.dump_observable(&BaseModel<VectorX>::get_energy_squared, "energy_squared", out_file);
    Gen.analyze_dataset("energy_squared", out_file, 100, 200, 200);
}

void DoMultiLevelMeasurementsFromFile(std::string filename) {

    HighFive::File file(std::string(DATA_DIR).append(filename), HighFive::File::ReadOnly);
    auto helper = file.getGroup("level0");//todo see if this step can be removed to be needed
    IsingModel test(helper);

    std::default_random_engine myengine{42L};
    MultiLevelHMCGenerator mygen(test, file, myengine);

    filename.insert(filename.rfind('/') + 1, "out_");
    std::string out_filename{std::string(DATA_DIR).append(filename)};
    DoMultiLevelMeasurements(mygen, out_filename);


}

void MultiLevelCriticalSimulation() {
    const int grid_size = 16;
    const int dim = 2;
    const int lambda = int_pow(grid_size, dim);
    const double C{0.1};
    const double beta{0.440686793509772};

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


    MultiLevelHMCGenerator mygen(test, {0, 16, 16}, {1, 16, 16}, 2, InterpolationType::Black_White, {8, 32, 64},
                                 {1. / 8., 1. / 32., 1. / 64.}, myengine);
    std::vector<double> acceptance_rates = mygen.generate_ensembles(phi0, 100000, 10000);
    for (auto acceptance_rate: acceptance_rates) {
        std::cout << "Acceptance rate:" << acceptance_rate << std::endl;
    }
    std::string filename{std::string(DATA_DIR).append(my_time).append(std::to_string(beta)).append(".h5")};
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    mygen.dumpToH5(file);

    std::string out_filename{std::string(DATA_DIR).append("out_").append(my_time).
            append(std::to_string(beta)).append(".h5")};
    DoMultiLevelMeasurements(mygen, out_filename);
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
    //test_hmc_measurements();
    //DoMultiLevelMeasurementsFromFile(
    // std::string("checker_board_multi_level_hmc_2_levels/02_05_2022__14_00_34_0.440687.h5"));
    MultiLevelCriticalSimulation();
    //return test_hip();
}


