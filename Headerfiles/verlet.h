//
// Created by cm on 06.05.21.
//

#ifndef __VERLET_H_
#define __VERLET_H_

void verletStep1(double &x, double &y, double &z, double &vx, double &vy, double &vz, double fx, double fy, double fz, double timestep);

void verletStep2(double &vx, double &vy, double &vz, double fx, double fy, double fz, double timestep);

void verletIntegratorConstantForce(double &xNow, double &yNow, double &zNow,
                                   double &vxNow, double &vyNow, double &vzNow,
                                   unsigned int nbSteps, double fx, double fy, double fz, double timestep);

#endif //__VERLET_H_
