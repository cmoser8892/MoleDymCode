//
// Created by cm on 13.05.21.
//

#ifndef MYPROJECT_ATOMS_H
#define MYPROJECT_ATOMS_H

#include "../Headerfiles/types.h"

class Atoms {
public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Names_t names;

    //set methods
    Atoms(const int &size) :
        positions{3,size}, velocities{3,size}, forces{3, size}
    {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
    }

    Atoms(const Positions_t &p) :
        positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}
    {
        velocities.setZero();
        forces.setZero();
    }

    Atoms(const Positions_t &p, const Velocities_t &v) :
        positions{p}, velocities{v}, forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
    }

    Atoms(const Names_t &n, const Positions_t &p, const Velocities_t &v) :
            positions{p}, velocities{v}, forces{3, p.cols()}, names{n} {
        assert(p.cols() == v.cols());
        forces.setZero();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }

};

#endif //MYPROJECT_ATOMS_H
