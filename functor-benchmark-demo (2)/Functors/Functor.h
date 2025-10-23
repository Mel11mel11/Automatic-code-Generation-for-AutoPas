#pragma once

template <typename Particle_T>
class Functor {
public:
    Functor() = default;

    virtual void AoSFunctor(Particle_T& p1, Particle_T& p2) = 0;
};