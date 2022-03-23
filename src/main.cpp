/**
  * @file       main.cpp
  * @brief      Starting point of the program
  * @author     nico
  * @version    0.0.1
  * @date       22.03.22
  */
#include <cmath>
#include <iostream>
#include <string>
#include <SrcConfig.h>
#include <Tensors_basic.h>


int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << argv[0] << " Version " << Tutorial_VERSION_MAJOR << "." << Tutorial_VERSION_MINOR << std::endl;
        std::cout << "Usage: " << argv[0] << " number" << std::endl;
        return 1;
    }
    // std::cout << SOURCE_DIR;
    // convert input to double
    const double inputValue = std::stod(argv[1]);

    // calculate square root

    const double outputValue = sqrt(inputValue);

    std::cout << "The square root of " << inputValue << " is " << outputValue
              << std::endl;
    return 0;
}
