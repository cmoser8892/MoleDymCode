//
// Created by cm on 06.05.21.
//
#include <iostream>
#include "../Headerfiles/verlet.h"

/**
 * @fn void verletStep1(double &x, double &y, double &z, double &vx, double &vy, double &vz, double fx, double fy, double fz, double timestep)
 * @brief implementation of the velocity-verlet algorithm
 * @param x
 * @param y
 * @param z
 * @param vx
 * @param vy
 * @param vz
 * @param fx
 * @param fy
 * @param fz
 * @param timestep
 */
void verletStep1(double & x, double & y, double & z, double & vx, double & vy, double & vz, double fx, double fy, double fz, double timestep) {
    //mass??? asume 1
    double mass = 1.0;
    //TODO pointer?? dont i have do derefence them
    //velocity
    vx += 0.5 * fx * timestep/mass;
    vy += 0.5 * fy * timestep/mass;
    vz += 0.5 * fz * timestep/mass;
    //position
    x += vx * timestep;
    y += vy * timestep;
    z += vz* timestep;
}

/**
 * @fn void verletStep2(double &vx, double &vy, double &vz, double fx, double fy, double fz, double timestep)
 * @brief implementation of the velocity-verlet algorithm
 * @param vx
 * @param vy
 * @param vz
 * @param fx
 * @param fy
 * @param fz
 * @param timestep
 */
void verletStep2(double &vx, double &vy, double &vz, double fx, double fy, double fz, double timestep) {
    //mass?? asume 1
    double mass = 1.0;
    vx += 0.5 * fx * timestep/mass;
    vy += 0.5 * fy * timestep/mass;
    vz += 0.5 * fz * timestep/mass;
}

/**
 * @fn verletIntegration over various steps
 * @brief Integration over varius steps
 * @param nbSteps
 * @param fx
 * @param fy
 * @param fz
 * @param timestep
 */
void verletIntegratorConstantForce(double &xNow, double &yNow, double &zNow,
                                   double &vxNow, double &vyNow, double &vzNow,
                                   unsigned int nbSteps, double fx, double fy, double fz, double timestep) {
    //
    double mass = 1; //gone anyway
    /** For loop */
    for(int i = 0; i < nbSteps; ++i) {
        //std::cout << "Step: " << i << std::endl;
        verletStep1(xNow,yNow,zNow,vxNow,vyNow,vzNow,fx,fy,fz,timestep);
        //TODO: compute force; constant anyway no to do anything need right?
        verletStep2(vxNow,vyNow,vzNow,fx,fy,fz,timestep);
    }
}