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
    //vars
    double totalPotentialEnergy = 0.0;
    int numberOfAtoms = atoms.nb_atoms();
    double delta = 0.0001;
    double thisAtomsPotential = 0.0;
    //TODO can i just discard all the past forces ?!
    atoms.forces = 0;
    //Calculate potential energy
    //TODO bad loop ?
    for(int i = 0; i < numberOfAtoms; ++i)
    {
        for(int j = i+1; j < numberOfAtoms; ++j)
        {
            /** Distance */
            //calculate the current energy between the atom i and j
            Vector_t vectorToOtherAtom = atoms.positions.col(j)-atoms.positions.col(i);
            double currentDistance = calculateDistanceBetweenVektors(vectorToOtherAtom);

            /** Energy */
            //calculate Lenard jones Potential
            thisAtomsPotential = calculateEnergy(currentDistance,sigma,epsilon);
            //add it up
            totalPotentialEnergy += thisAtomsPotential;

            /** Forces */
            //normalize vector
            Vector_t normalizedVectorToOtherAtom = vectorToOtherAtom/currentDistance;
            //calculate the deltaV
            //TODO Taylor??
            double deltaThisAtomsPotential= calculateEnergy(currentDistance+delta,epsilon,sigma) - calculateEnergy(currentDistance - delta, epsilon,sigma);
            //put the force there
            atoms.forces.col(i) += 68* (deltaThisAtomsPotential/(2*delta)) * normalizedVectorToOtherAtom;
            //other atom has the force in the other direction
            atoms.forces.col(j) += -1* atoms.forces.col(i);
        }
    }
    //first part of equation TODO do i even need that as the energy is not counted twice?!
    //totalPotentialEnergy*= 0.5;
    //std::cout << atoms.forces << std::endl;
    return thisAtomsPotential;
}

/**
 * @fn double calculateDistanceBetweenVektors(Vector_t distanceVector)
 * @brief calculates the length of the Vector
 * @param distanceVector
 * @return double the length of the Vektor
 */
double calculateDistanceBetweenVektors(Vector_t distanceVector) {
    //pythagoras
    return sqrt(distanceVector(0)* distanceVector(0) +
                                   distanceVector(1)* distanceVector(1) +
                                   distanceVector(2)* distanceVector(2));
}

/**
 * @fn double calculateEnergy(double distance, double epsilon, double sigma)
 * @brief calculates the Potential Energy between two atoms V(r) using the Lenard Jones Potential
 * @param distance
 * @param epsilon
 * @param sigma
 * @return double the Potential energy
 */
double calculateEnergy(double distance, double epsilon, double sigma)
{
    return 4*epsilon*(pow(sigma/distance,12) - pow(sigma/distance,6));
}