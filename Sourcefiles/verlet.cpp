//
// Created by cm on 06.05.21.
//
#include <iostream>
#include "../Headerfiles/verlet.h"
#include "../Headerfiles/types.h"

void verletStep1(double & x, double & y, double & z, double & vx, double & vy, double & vz, double fx, double fy, double fz, double timestep) {
    //mass??? asume 1
    double mass = 1.0;
    //velocity
    vx += 0.5 * fx * timestep/mass;
    vy += 0.5 * fy * timestep/mass;
    vz += 0.5 * fz * timestep/mass;
    //position
    x += vx * timestep;
    y += vy * timestep;
    z += vz* timestep;
}

void verletStep1(Positions_t &positions, Velocities_t &velocities, const Forces_t &forces ,double timestep)
{
    double mass = 1.0;
    //velocity
    velocities += 0.5 * forces * timestep/mass;
    //postions
    positions += velocities * timestep;
}


void verletStep2(double &vx, double &vy, double &vz, double fx, double fy, double fz, double timestep) {
    //mass?? asume 1
    double mass = 1.0;
    vx += 0.5 * fx * timestep/mass;
    vy += 0.5 * fy * timestep/mass;
    vz += 0.5 * fz * timestep/mass;
}

void verletStep2(Velocities_t &velocities, Forces_t &forces, double timestep)
{
    double mass = 1.0;
    velocities += 0.5 * forces * timestep/mass;
}


void verletIntegratorConstantForce(double &xNow, double &yNow, double &zNow,
                                   double &vxNow, double &vyNow, double &vzNow,
                                   unsigned int nbSteps, double fx, double fy, double fz, double timestep) {
    //
    double mass = 1; //gone anyway
    /** For loop */
    for(int i = 0; i < nbSteps; ++i) {
        //std::cout << "Step: " << i << std::endl;
        //verletStep1(xNow,yNow,zNow,vxNow,vyNow,vzNow,fx,fy,fz,timestep);
        //TODO: compute force; constant anyway no to do anything need right?
        //verletStep2(vxNow,vyNow,vzNow,fx,fy,fz,timestep);
    }
}

void verletIntegratorConstantForce(Positions_t &positions, Velocities_t &velocities, Forces_t &forces
        , double timestep, unsigned int nbSteps)
{
    double mass = 1;
    // For loop
    for(int i = 0; i < nbSteps; ++i) {
        //std::cout << "Step: " << i << std::endl;
        verletStep1(positions,velocities,forces,timestep);
        verletStep2(velocities,forces,timestep);
    }
}