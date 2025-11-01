#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include "../Particle.h"
#include "../Functors/LJFunctorGenerated.h"
#include "../Functors/LJFunctorReference.h"
#include "../generated_files/generated_force.hpp"
#include "../Functors/GravFunctorReference.h"
#include "../Functors/GravFunctorGenerated.h"
#include "../Functors/MieFunctorGenerated.h"
#include "../Functors/MieFunctorReference.h"

std::vector<Particle> makeGrid(int N, double spacing) {
    std::vector<Particle> ps;
    ps.reserve(N * N * N);
    int id = 0;
    for (int x = 0; x < N; ++x)
        for (int y = 0; y < N; ++y)
        for (int z = 0; z < N; ++z) {
            Particle p;
            p.setID(id++);
            p.setR({x * spacing, y * spacing, z * spacing});
            p.setF({0, 0, 0});
            ps.push_back(p);
        }
    return ps;
}

template <class F, class P>
void runPairs(std::vector<P>& ps, F& functor) { // function from main.cpp
    const size_t N = ps.size();
    for (size_t i = 0; i < N; ++i)
        for (size_t j = i + 1; j < N; ++j)
            functor.AoSFunctor(ps[i], ps[j]);
}

std::array<double, 3> sumForces(const std::vector<Particle>& ps) { 
    std::array<double, 3> s{0, 0, 0};
    for (auto& p : ps) {
        const auto& f = p.getF();
        s[0] += f[0]; s[1] += f[1]; s[2] += f[2];
    }
    return s;
}
auto copy(std::vector<Particle> grid){
  auto neu_grid = grid;
return neu_grid;

}
// space for main test
int main() {
    int N = 4; // 4x4x4 grid
    double distance = 1.2; // On the advice of Markus
    double sig = 1.0, eps = 1.0;
    const int n=7, m=6;
    const double G=6.67430e-11;
    bool newton3 = true;
    auto ps_ref = makeGrid(N, distance);
    auto ps_gen = copy(ps_ref); // use the same grid
    // create the forces
    LJFunctorReference<Particle> ref_lj(sig, eps, newton3);
    LJFunctorGenerated<Particle> gen_lj(sig, eps, newton3);

    MieFunctorGenerated<Particle> miegen(sig, eps, n, m,  newton3);
    MieFunctorReference<Particle> mieRef(sig, eps, n, m,  newton3);
    GravFunctorGenerated<Particle> gen_grav(G, newton3);
    GravFunctorReference<Particle> ref_grav(G, newton3);
    // calculate the forces
    runPairs(ps_ref, ref_lj);
    runPairs(ps_gen, gen_lj);
    runPairs(ps_ref, mieRef);
    runPairs(ps_gen, miegen);
    runPairs(ps_ref, ref_grav);
    runPairs(ps_gen, gen_grav);
auto ps_lj_ref  = ps_ref;
auto ps_lj_gen  = ps_gen;
auto ps_gr_ref  = ps_ref;
auto ps_gr_gen  = ps_gen;
auto ps_mie_ref = ps_ref;
auto ps_mie_gen = ps_gen;

    // calculate the differences
    /*double maxDelta = 0.0, l2 = 0.0;
    for (size_t i = 0; i < ps_ref.size(); ++i) { // ps_ref.size() tells us how many particles we have 4x4x4
        auto f1 = ps_ref[i].getF();
        auto f2 = ps_gen[i].getF();
        for (int k = 0; k < 3; ++k) {
            double d = f1[k] - f2[k];
            maxDelta = std::max(maxDelta, std::abs(d)); // |d| 
            //l2 += d * d;
        }
    }*/
   //Lennard-Jones 
    double maxDelta_LJ = 0.0;
    for (size_t i=0;i<ps_lj_ref.size();++i)
    for (int k=0;k<3;++k)
    maxDelta_LJ = std::max(maxDelta_LJ, std::abs(ps_lj_ref[i].getF()[k]-ps_lj_gen[i].getF()[k]));

   //Gravity 
    double maxDelta_Grav = 0.0;
    for (size_t i=0;i<ps_gr_ref.size();++i)
    for (int k=0;k<3;++k)
    maxDelta_Grav = std::max(maxDelta_Grav, std::abs(ps_gr_ref[i].getF()[k]-ps_gr_gen[i].getF()[k])); 
   //l2 = std::sqrt(l2 / ps_ref.size());
   //Mie
    double maxDelta_Mie = 0.0;  
    for (size_t i=0;i<ps_mie_ref.size();++i)
    for (int k=0;k<3;++k)
    maxDelta_Mie = std::max(maxDelta_Mie, std::abs(ps_mie_ref[i].getF()[k]-ps_mie_gen[i].getF()[k]));

runPairs(ps_lj_ref,  ref_lj);
runPairs(ps_lj_gen,  gen_lj);
runPairs(ps_gr_ref,  ref_grav);
runPairs(ps_gr_gen,  gen_grav);
runPairs(ps_mie_ref, mieRef);
runPairs(ps_mie_gen, miegen);

// fark ve toplamları potansiyel bazında hesapla
auto Lj_sRef  = sumForces(ps_lj_ref);
auto Lj_sGen  = sumForces(ps_lj_gen);
//auto Grav_sRef= sumForces(ps_gr_ref);
//auto Grav_sGen= sumForces(ps_gr_gen);
// auto Mie_sRef = sumForces(ps_mie_ref);
// auto Mie_sGen = sumForces(ps_mie_gen);

    std::cout << "Newton3=" << newton3 << "\n";
    std::cout << "Lennard-Jones ΣF(Ref) = (" << Lj_sRef[0] << ", " << Lj_sRef[1] << ", " << Lj_sRef[2] << ")\n";
    std::cout << "Lennard-Jones ΣF(Gen) = (" << Lj_sGen[0] << ", " << Lj_sGen[1] << ", " << Lj_sGen[2] << ")\n";
    std::cout << "max|ΔF| = " << maxDelta_LJ << "\n";
    //std:: cout<<  "mean L2 = " << l2 << "\n";
    // std::cout << "Grav ΣF(Ref) = (" << Grav_sRef[0] << ", " << Grav_sRef[1] << ", " << Grav_sRef[2] << ")\n";
   // std::cout << "Grav ΣF(Gen) = (" << Grav_sGen[0] << ", " << Grav_sGen[1] << ", " << Grav_sGen[2] << ")\n";
    // std::cout << "max|ΔF| = " << maxDelta_Grav << "\n";
    //std:: cout<<  "mean L2 = " << l2 << "\n";
   // std::cout << "Mie ΣF(Ref) = (" << Mie_sRef[0] << ", " << Mie_sRef[1] << ", " << Mie_sRef[2] << ")\n";
    // std::cout << "Mie ΣF(Gen) = (" << Mie_sGen[0] << ", " << Mie_sGen[1] << ", " << Mie_sGen[2] << ")\n";
    // std::cout << "max|ΔF| = " << maxDelta_Mie << "\n";
    //std:: cout<<  "mean L2 = " << l2 << "\n";
    return 0;
}
