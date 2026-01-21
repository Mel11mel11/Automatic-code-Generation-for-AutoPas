#pragma once
#include <array>
#include "ArrayMath.h"

class Particle {
public:
    Particle() = default;
double getFx() const { return _force[0]; }
double getFy() const { return _force[1]; }
double getFz() const { return _force[2]; }

    // position, velocity, force, id, (opsiyonel) mass
    Particle(std::array<double,3> r,
             std::array<double,3> v,
             std::array<double,3> f,
             int id,
             double mass = 1.0)
        : _position(r), _velocity(v), _force(f), _id(id), _mass(mass) {}

    // ---- Force accessors ----
    const std::array<double,3>& getF() const { return _force; }
    void setF(std::array<double,3> f) { _force = f; }

    // Particle.h  
    void addF(std::array<double,3> fAdd) {
    _force[0] += fAdd[0];
    _force[1] += fAdd[1];
    _force[2] += fAdd[2];
        }
    void subF(std::array<double,3> fSub) {
    _force[0] -= fSub[0];
    _force[1] -= fSub[1];
    _force[2] -= fSub[2];
    }


    const std::array<double,3>& getR() const { return _position; }

    void setR(std::array<double,3> r) { _position = r; }

    void setID(size_t id) { _id = id; }

    double getMass() const { return _mass; }

    void setMass(double m) { _mass = m; }

private:
    std::array<double,3> _position{};
    std::array<double,3> _velocity{};
    std::array<double,3> _force{};
    size_t _id{};
    double _mass{1.0};   // 
};
