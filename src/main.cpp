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
#include <LeapFrogIntegrator.h>

void test_leap_frog() {
    const int grid_size = 8;
    const int dim = 2;
    const int lambda = int_pow(grid_size, dim);
    const double C{4.1};
    const double beta{2.0};

    VectorX phi0(lambda);
    phi0.setRandom();
    VectorX pi0(lambda);
    pi0.setRandom();
    VectorX h0(lambda);
    h0.setRandom();
    VectorX eta0(lambda);
    eta0.setZero();

    IsingModel test(beta, h0, eta0, C, dim, 1, grid_size);

    LeapFrogIntegrator leapTest(test);
    int numMD = 10;

    double S_start = test.get_action(phi0) + pi0.dot(pi0) / 2.;
    std::vector<double> S_error;
    std::vector<int> MD;
    for (int i = 10; i < 200; i += 20) {
        MD.push_back(i);
        VectorX phiNew(phi0);
        VectorX piNew(pi0);
        std::cout << 1. / i << std::endl;
        leapTest.integrate(i, 1. / i, phiNew, piNew);
        double S_end = test.get_action(phiNew) + piNew.dot(piNew) / 2.;
        S_error.push_back(abs(S_start - S_end) / S_start);
        std::cout << i << " : " << abs(S_start - S_end) / S_start << std::endl;
    }

}


int main(int argc, char *argv[]) {
    std::cout << "Hello\n";
    test_leap_frog();
}
