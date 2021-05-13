//
// Created by cm on 13.05.21.
//

#include "../Headerfiles/lenardJonesDirectionSummation.h"
#include <iostream>
/**
 * @func double lendardJonesDirectSummation(Atoms &atoms, double epsilon, double sigma)
 * @brief calculates the Lenard Jones Potential as it is specificated in Milestone 4
 * @param atoms Datastruktur of all the atoms in the system; directly modifies the force
 * @param epsilon Part of equation
 * @param sigma Part of equation
 * @return the Potential energy
 */
double lendardJonesDirectSummation(Atoms &atoms, double epsilon, double sigma) {
    double totalPotentialEnergy = 0.0;
    int numberOfAtoms = atoms.nb_atoms();
    //TODO can i just discard all the past forces ?!
    atoms.forces = 0;
    //Calculate potential energy
    for(int i = 0; i < numberOfAtoms; ++i)
    {
        for(int j = i+1; j < numberOfAtoms; ++j)
        {
            //calculate the current energy in the system
            Vector_t pointerToOtherAtom = atoms.positions.col(j)-atoms.positions.col(i);
            //std::cout << pointerToOtherAtom << std::endl;
            double currentDistance = calculateDistanceBetweenVektors(pointerToOtherAtom);
            double thisAtomsPotential = 4*epsilon*(pow(sigma/currentDistance,12) + pow(sigma/currentDistance,6));
            //std::cout << thisAtomsPotential<< std::endl;
            totalPotentialEnergy += thisAtomsPotential;
            //forces
            Vector_t normalizedPointerToOtherAtom = pointerToOtherAtom/currentDistance;
            //std::cout << normalizedPointerToOtherAtom  << std::endl;
            //TODO: check here!
            atoms.forces.col(i) += totalPotentialEnergy/currentDistance * normalizedPointerToOtherAtom;
            //std::cout << atoms.forces.col(i)  << std::endl;
            //other atom has the force in the other direction
            atoms.forces.col(j) += -1* atoms.forces.col(i);
        }
    }
    //first part of equation TODO do i even need that as the energy is not counted twice?!
    totalPotentialEnergy*= 0.5;
    //std::cout << atoms.forces << std::endl;
    return totalPotentialEnergy;
}

double calculateDistanceBetweenVektors(Vector_t distanceVector) {
    return sqrt(distanceVector(0)* distanceVector(0) +
                                   distanceVector(1)* distanceVector(1) +
                                   distanceVector(2)* distanceVector(2));
}