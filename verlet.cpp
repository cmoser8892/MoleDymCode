//
// Created by cm on 06.05.21.
//
#include "verlet.h"

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
void verletStep1(double & x, double & y, double & z, double & vx, double & vy, double & vz, double fx, double fy, double fz, double timestep){
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
void verletStep2(double &vx, double &vy, double &vz, double fx, double fy, double fz, double timestep){
    //mass?? asume 1
    double mass = 1.0;
    vx += 0.5 * fx * timestep/mass;
    vy += 0.5 * fy * timestep/mass;
    vz += 0.5 * fz * timestep/mass;
}