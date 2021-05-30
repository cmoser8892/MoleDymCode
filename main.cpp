#include <iostream>

#include "Headerfiles/verlet.h"
#include "Headerfiles/types.h"
#include "Headerfiles//lenardJonesDirectionSummation.h"
#include "Headerfiles/xyz.h"
#include<fstream>

using namespace std;

int main() {
    std::cout << "Molecular Dynamics Project" << std::endl;
    auto [names, positions, velocities]{read_xyz_with_velocities("../AJupyter/lj54.xyz")};
    return 0;
}
