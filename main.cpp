#include <iostream>

#include "verlet.h"

int main() {
    std::cout << "Hello world!" << std::endl;
    std::cout << "TEST" << std::endl;
    double testNumber = 1.0;
    verletStep1(testNumber, testNumber, testNumber, testNumber, testNumber, testNumber, testNumber, testNumber, testNumber, testNumber);
    verletStep2(testNumber, testNumber, testNumber, testNumber, testNumber, testNumber, testNumber);
    std::cout << testNumber << std::endl;
    return 0;
}
