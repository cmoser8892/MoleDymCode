//
// Created by cm on 06.05.21.
//

#ifndef __VERLET_H
#define __VERLET_H

//TODO test for todo big nice

void verletStep1(double &x, double &y, double &z, double &vx, double &vy, double &vz, double fx, double fy, double fz, double timestep);

void verletStep2(double &vx, double &vy, double &vz, double fx, double fy, double fz, double timestep);

#endif //__VERLET_H
