#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

#include "Particle.h"

// Functors
#include "../Functors/LJFunctorReference.h"
#include "../Functors/generatednew /generated_LJ_O000.hpp"

#include "../Functors/MieFunctorReference.h"
#include "../Functors/generatednew /generated_Mie_O000.hpp"

#include "../Functors/GravFunctorReference.h"
#include "../Functors/generatednew /generated_Gravity_O000.hpp"

#include "../Functors/KryptonFunctorReference.h"
#include "../Functors/generatednew /generated_Krypton_O000.hpp"


std::vector<Particle> makeGrid(int N, double spacing) {
    std::vector<Particle> ps;
    ps.reserve(N*N*N);
    int id = 0;

    for(int x=0;x<N;x++)
        for(int y=0;y<N;y++)
            for(int z=0;z<N;z++) {
                Particle p;
                p.setID(id++);
                p.setR({x*spacing, y*spacing, z*spacing});
                p.setF({0,0,0});
                ps.push_back(p);
            }

    return ps;
}


template<class F,class P>
void runPairs(std::vector<P>& ps, F& functor) {
    size_t N = ps.size();
    for(size_t i=0;i<N;i++)
        for(size_t j=i+1;j<N;j++)
            functor.AoSFunctor(ps[i], ps[j]);
}

std::array<double,3> sumForces(const std::vector<Particle>& ps) {
    std::array<double,3> s{0,0,0};
    for(auto& p: ps) {
        auto f = p.getF();
        s[0]+=f[0]; s[1]+=f[1]; s[2]+=f[2];
    }
    return s;
}

template<class P>
double maxDelta(const std::vector<P>& R, const std::vector<P>& G) {
    double m=0.0;
    for(size_t i=0;i<R.size();i++) {
        auto f1 = R[i].getF();
        auto f2 = G[i].getF();
        for(int k=0;k<3;k++)
            m = std::max(m, std::abs(f1[k]-f2[k]));
    }
    return m;
}

template<class REF,class GEN>
void testPotential(const std::string& name,
                   REF& refF, GEN& genF,
                   int N, double spacing)
{
    auto ps_ref = makeGrid(N, spacing);
    auto ps_gen = ps_ref;

    runPairs(ps_ref, refF);
    runPairs(ps_gen, genF);

    auto Sref = sumForces(ps_ref);
    auto Sgen = sumForces(ps_gen);
    double dmax = maxDelta(ps_ref, ps_gen);

    std::cout << "\n=== " << name << " ===\n";
    std::cout << "ΣF(Ref) = (" << Sref[0] << ", " << Sref[1] << ", " << Sref[2] << ")\n";
    std::cout << "ΣF(Gen) = (" << Sgen[0] << ", " << Sgen[1] << ", " << Sgen[2] << ")\n";
    std::cout << "max|ΔF| = " << dmax << "\n";
}

int main() {

    const int N = 4;
    const double spacing = 1.2;
    const bool newton3 = false;

    const double sigma = 1.0;
    const double epsilon = 1.0;
    const int n = 7, m = 6;
    const double G = 6.67430e-11;

    // Krypton params --------------------
    const double A   = 5.2239e3;
    const double a1  = 3.2909;
    const double a2  = -3.1237e-1;
    const double a_m1 = 5.6840;
    const double b    = 3.0;
    const double C6   = 2.451e2;
    const double C8   = 1.684e4;
    const double C10  = 1.618e6;
    const double cutoff = 2.5;
    std::cout << "Newton3 = " << newton3 << "\n";

    // Lennard-Jones
    LJFunctorReference<Particle> ljR(sigma,epsilon,newton3,cutoff);
    LJFunctor_Gen_O000<Particle> ljG(sigma,epsilon,newton3,cutoff);
    
    testPotential("Lennard-Jones", ljR, ljG, N, spacing);

    // Mie
    MieFunctorReference<Particle> mieR(sigma,epsilon,n,m,newton3);
    MieFunctor_Gen_O000<Particle> mieG(sigma,epsilon,n,m,newton3);
    testPotential("Mie Potential", mieR, mieG, N, spacing);

    // Gravity
    GravFunctorReference<Particle> grR(G,newton3);
    GravityFunctor_Gen_O000<Particle> grG(G,newton3);
    testPotential("Gravity", grR, grG, N, spacing);

    // Krypton
    KryptonFunctorReference<Particle> kR(A,a1,a2,a_m1,b,C6,C8,C10,newton3);
    KryptonFunctorGenerated_Gen_O000<Particle> kG(A,a1,a2,a_m1,b,C6,C8,C10,newton3);
    testPotential("Krypton", kR, kG, N, spacing);

    return 0;
}
