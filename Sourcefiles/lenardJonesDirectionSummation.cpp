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
            double currentDistance = calculateDistanceBetweenVektors(vectorToOtherAtom);

            /** Energy */
            totalPotentialEnergy += calculateEnergy(currentDistance,epsilon,sigma);

            /** Forces */
            //normalize vector
            Vector_t normalizedVectorToOtherAtom = vectorToOtherAtom/currentDistance;
            //calculate the deltaV
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
    return totalPotentialEnergy ;
}

/**
 * @fn double lenardJonesDirectSummationWithCutoff(Atoms &atoms, NeighborList &list, double epsilon, double sigma)
 * @brief updated lenardJonesDirectSummation updated for Milestone 6
 * @param atoms
 * @param list
 * @param epsilon
 * @param sigma
 * @return potential Energy in the system
 */
double lenardJonesDirectSummationWithCutoff(Atoms &atoms, NeighborList &list, double epsilon, double sigma) {
    //basic vars
    double totalPotentialEnergy = 0.0;
    atoms.forces.setZero();
    list.update(atoms);
    //loop
    for(auto[i,j]:list) {
        if(i <j) {
            /** Calculate a the distance vector*/
            Vector_t vectorToOtherAtom = atoms.positions.col(j)-atoms.positions.col(i);
            double currentDistance = calculateDistanceBetweenVektors(vectorToOtherAtom);

            /** Energy calculation */
            totalPotentialEnergy += calculateEnergy(currentDistance,epsilon,sigma);

            /** Force calculation */
            //normalize vector
            Vector_t normalizedVectorToOtherAtom = vectorToOtherAtom/currentDistance;
            //calculate the deltaV
            double deltaForce = calculateForceAnalytical(epsilon,sigma,currentDistance);
            Vector_t force = deltaForce * normalizedVectorToOtherAtom;
            //put the force there
            atoms.forces.col(i) += force;
            //other atom has the force in the other direction
            atoms.forces.col(j) -= force;
        }
    }
    //
    return totalPotentialEnergy;
}

/** Satans little Helpers */
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
    return 4.0*epsilon*(pow(sigma/distance,12.) - pow(sigma/distance,6.));
}

/**
 * @fn double calculateForceAnalytical(double epsilon, double sigma, double distance)
 * @brief calculates the force analytically form the derivation of the Lenard Jones Potential
 * @param epsilon
 * @param sigma
 * @param distance
 * @return force
 */
double calculateForceAnalytical(double epsilon, double sigma, double distance)
{
    return 4* epsilon *((6*pow(sigma,6)/pow(distance,7)) - (12*pow(sigma,12)/pow(distance,13)));
}


