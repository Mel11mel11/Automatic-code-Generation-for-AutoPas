#pragma once
#include <array>
#include "ArrayMath.h"

class Particle {
public:
    Particle() = default;

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

    void addF(std::array<double,3> fAdd) {
        using namespace arrayMath::literals;
        _force += fAdd;
    }
    void subF(std::array<double,3> fSub) {
        using namespace arrayMath::literals;
        _force -= fSub;
    }

    // ---- Position accessors ----
    const std::array<double,3>& getR() const { return _position; }
    void setR(std::array<double,3> r) { _position = r; }

    // ---- Mass accessors (YENİ) ----
    double getMass() const { return _mass; }

    void   setMass(double m) { _mass = m; }

private:
    std::array<double,3> _position{};
    std::array<double,3> _velocity{};
    std::array<double,3> _force{};
    size_t _id{};
    double _mass{1.0};   // <-- mutlaka ekli olmalı
};
