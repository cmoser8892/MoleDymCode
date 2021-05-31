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
    //forces back to 0/discard old forces
    atoms.forces.setZero();
    //Calculate potential energy
    for(int i = 0; i < numberOfAtoms; ++i)
    {
        for(int j = i+1; j < numberOfAtoms; ++j)
        {
            /** Distance */
            //calculate the current energy between the atom i and j
            Vector_t vectorToOtherAtom = atoms.positions.col(j)-atoms.positions.col(i);
            //std::cout << vectorToOtherAtom << std::endl;
            double currentDistance = calculateDistanceBetweenVektors(vectorToOtherAtom);

            /** Energy */
            totalPotentialEnergy += calculateEnergy(currentDistance,epsilon,sigma);

            /** Forces */
            //normalize vector
            Vector_t normalizedVectorToOtherAtom = vectorToOtherAtom/currentDistance;
            //calculate the deltaV
            //Vector_t deltaForce = calculateForceAnalytical(epsilon,sigma,atoms.positions.col(j),atoms.positions.col(i)); //placeholder
            double deltaForce = calculateForceAnalytical(epsilon,sigma,currentDistance);
            Vector_t force = deltaForce * normalizedVectorToOtherAtom;
            //put the force there
            atoms.forces.col(i) += force;
            //other atom has the force in the other direction
            atoms.forces.col(j) -= force;
        }
    }
    //first part of equation
    //totalPotentialEnergy*= 0.5;
    //std::cout << atoms.forces << std::endl;
    return totalPotentialEnergy ;
}

/**
 * @fn double calculateDistanceBetweenVektors(Vector_t distanceVector)
 * @brief calculates the length of the Vector
 * @param distanceVector
 * @return double the length of the Vektor
 */
double calculateDistanceBetweenVektors(Vector_t distanceVector) {
    //pythagoras
    double dist = distanceVector(0)* distanceVector(0) +
                  distanceVector(1)* distanceVector(1) +
                  distanceVector(2)* distanceVector(2);
    return sqrt(dist);
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
    double first = pow(sigma/distance,12.);
    double second = pow(sigma/distance,6.);
    return 4.0*epsilon*(first - second);
}

/**
 * @fn
 * @brief
 * @param epsilon
 * @param sigma
 * @param distance
 * @return
 */
double calculateForceAnalytical(double epsilon, double sigma, double distance)
{
    double first = 6*pow(sigma,6)/pow(distance,7);
    double second = 12*pow(sigma,12)/pow(distance,13);
    return 4* epsilon *(first - second);
}