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
    IsingModel test(3,10,2,0.);
    test.print_connectivity_matrix();
}
