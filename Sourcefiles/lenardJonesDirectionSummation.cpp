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
    //forces back to 0/discard old forces
    atoms.forces = 0;
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
            //calculate Lenard jones Potential
            thisAtomsPotential = calculateEnergy(currentDistance,sigma,epsilon);
            //add it up
            totalPotentialEnergy += thisAtomsPotential;

            /** Forces */
            //normalize vector
            Vector_t normalizedVectorToOtherAtom = vectorToOtherAtom/currentDistance;
            //calculate the deltaV
            double deltaForce = calculateForceAnalytical(epsilon,sigma,atoms.positions.col(i),atoms.positions.col(j)); //placeholder
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

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma){
    double E = 0.;
    double ssq = sigma*sigma;
    for(int j = 0; j<atoms.nb_atoms(); j++) {// Looping over the upper triangle of the pair matrix
        for (int i = j+1; i < atoms.nb_atoms(); i++) {
            Vector_t rij_v = atoms.positions.col(j)-atoms.positions.col(i);
            double rij_sq = rij_v(0)*rij_v(0)+rij_v(1)*rij_v(1)+rij_v(2)*rij_v(2);
            E += 4.*epsilon*(pow(ssq/rij_sq,6.)-pow(ssq/rij_sq,3.)); //no 1/2, because E = Eij+Eji
            Vector_t Fij = 24.*epsilon*(pow(ssq*rij_sq,3.)-2.*pow(ssq,6.))/pow(rij_sq,7.)*rij_v;
            atoms.forces.col(i) += Fij;
            atoms.forces.col(j) -= Fij;
        }
    }
    return E;
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
    return 4*epsilon*(pow(sigma/distance,12.) - pow(sigma/distance,6.));
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
double calculateForceAnalytical(double epsilon, double sigma, Vector_t vectorI, Vector_t vectorJ) {
    //TODO vector or double?
    double returnValue = 0;
    //do in steps to avoid confusion
    /** Differences in X,Y and Z for later  */
    double xMinusIJ = vectorI(0) - vectorJ(0);
    double yMinusIJ = vectorI(1) - vectorJ(1);
    double zMinusIJ = vectorI(2) - vectorJ(2);

    /** Denominators for the Calculations s*/
    double denominatorsInnerSum = pow(xMinusIJ,2)
                                  +pow(yMinusIJ,2)
                                  +pow(zMinusIJ,2);
    double denomniatorPart1 = pow(denominatorsInnerSum,7.);
    double denomniatorPart2 = pow(denominatorsInnerSum,4.);

    /** Nominator for the calc */
    //TODO delta ?? 2dik - 2djk first part always 1 second part always 0?!
    //TODO i only need to compute this for i = k right j = k is computed in the function already right?
    double nominatorPart1Common = 6*pow(sigma,12.) * 2;
    double nominatorPart2Common = 3*pow(sigma,6.)  * 2;

    /** putting the inner equations together */
    //this could be written more intelligent but for readability it stays like this; xyzMinusIJ placed more intelligent
    double innerEquationX = - ((nominatorPart1Common*xMinusIJ)/denomniatorPart1) + ((nominatorPart2Common*xMinusIJ)/denomniatorPart2);
    double innerEquationY = - ((nominatorPart1Common*yMinusIJ)/denomniatorPart1) + ((nominatorPart2Common*yMinusIJ)/denomniatorPart2);
    double innerEquationZ = - ((nominatorPart1Common*zMinusIJ)/denomniatorPart1) + ((nominatorPart2Common*zMinusIJ)/denomniatorPart2);

    /** putting the whole equation together */
    double innerEquation =  innerEquationX + innerEquationY + innerEquationZ;
    returnValue = 4*epsilon*  innerEquation;
    /** */
    return returnValue;
}