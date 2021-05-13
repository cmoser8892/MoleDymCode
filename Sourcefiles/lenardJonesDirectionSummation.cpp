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
    double potentialEnergy = 0.0;
    int numberOfAtoms = atoms.nb_atoms();
    //Calculate potential energy
    for(int i = 0; i < numberOfAtoms; ++i)
    {
        for(int j = i+1; j < numberOfAtoms; ++j)
        {
            double currentDistance = calculateDistanceBetweenVektors(atoms.positions.col(i)-atoms.positions.col(j));
            potentialEnergy += 4*epsilon*(
                    pow(sigma/currentDistance,12) + pow(sigma/currentDistance,6)
                    );
        }
    }
    //first part of equation
    potentialEnergy *= 0.5;

    return potentialEnergy;
}

double calculateDistanceBetweenVektors(Vector_t distanceVector) {
    return sqrt(distanceVector(0)* distanceVector(0) +
                                   distanceVector(1)* distanceVector(1) +
                                   distanceVector(2)* distanceVector(2));
}