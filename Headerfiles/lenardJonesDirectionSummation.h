//
// Created by cm on 13.05.21.
//

#ifndef MYPROJECT_LENARDJONESDIRECTIONSUMMATION_H
#define MYPROJECT_LENARDJONESDIRECTIONSUMMATION_H

#include "../Headerfiles/atoms.h"

double lendardJonesDirectSummation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

double calculateDistanceBetweenVektors(Vector_t distanceVector);
double calculateEnergy(double distance, double epsilon, double sigma);

double calculateForceAnalytical(double epsilon, double sigma, Vector_t vectorI, Vector_t vectorJ);


#endif //MYPROJECT_LENARDJONESDIRECTIONSUMMATION_H
