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
            double deltaForce = calculateForceAnaly(epsilon,sigma,currentDistance);
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
 * @param vectorI
 * @param vectorJ
 * @return
 */
Vector_t calculateForceAnalytical(double epsilon, double sigma, Vector_t vectorI, Vector_t vectorJ) {
    Vector_t returnValue;
    //do in steps to avoid confusion
    /** Differences in X,Y and Z for later  */
    double xMinusIJ = vectorI(0) - vectorJ(0);
    double yMinusIJ = vectorI(1) - vectorJ(1);
    double zMinusIJ = vectorI(2) - vectorJ(2);

    /** Denominators for the Calculations s*/
    double denominatorsInnerSum =  pow(xMinusIJ,2)
                                  +pow(yMinusIJ,2)
                                  +pow(zMinusIJ,2);
    double denomniatorPart1 = pow(denominatorsInnerSum,7);
    double denomniatorPart2 = pow(denominatorsInnerSum,4);

    /** Nominator for the calc */
    double nominatorPart1Common = 6*pow(sigma,12)   *2;
    double nominatorPart2Common = 3*pow(sigma,6)    *2;

    /** putting the inner equations together */
    double innerEquationX = - ((nominatorPart1Common*xMinusIJ)/denomniatorPart1) + ((nominatorPart2Common*xMinusIJ)/denomniatorPart2);
    double innerEquationY = - ((nominatorPart1Common*yMinusIJ)/denomniatorPart1) + ((nominatorPart2Common*yMinusIJ)/denomniatorPart2);
    double innerEquationZ = - ((nominatorPart1Common*zMinusIJ)/denomniatorPart1) + ((nominatorPart2Common*zMinusIJ)/denomniatorPart2);

    /** putting the whole equation together */
    Vector_t innerEquation(innerEquationX , innerEquationY ,innerEquationZ);
    returnValue = 4 * epsilon* innerEquation;
    /** */
    return returnValue;
}

double calculateForceAnaly(double epsilon, double sigma, double distance)
{
    double first = 6*pow(sigma,6)/pow(distance,7);
    double second = 12*pow(sigma,12)/pow(distance,13);
    return 4* epsilon *(first - second);
}