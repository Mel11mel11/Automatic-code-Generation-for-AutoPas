#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include "../Particle.h"
#include "../Functors/LJFunctorGenerated.h"
#include "../Functors/LJFunctorReference.h"
#include "../generated_files/generated_force.hpp"

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
    bool newton3 = true;
    auto ps_ref = makeGrid(N, distance);
    auto ps_gen = copy(ps_ref); // use the same grid
    // create the forces
    LJFunctorReference<Particle> ref_lj(sig, eps, newton3);
    LJFunctorGenerated<Particle> gen_lj(sig, eps, newton3);
    // calculate the forces
    runPairs(ps_ref, ref_lj);
    runPairs(ps_gen, gen_lj);
    // calculate the differences
    double maxDelta = 0.0, l2 = 0.0;
    for (size_t i = 0; i < ps_ref.size(); ++i) { // ps_ref.size() tells us how many particles we have 4x4x4
        auto f1 = ps_ref[i].getF();
        auto f2 = ps_gen[i].getF();
        for (int k = 0; k < 3; ++k) {
            double d = f1[k] - f2[k];
            maxDelta = std::max(maxDelta, std::abs(d)); // |d| 
            //l2 += d * d;
        }
    }
   // l2 = std::sqrt(l2 / ps_ref.size());    
    auto sRef = sumForces(ps_ref);
    auto sGen = sumForces(ps_gen);
    std::cout << "Newton3=" << newton3 << "\n";
    std::cout << "ΣF(Ref) = (" << sRef[0] << ", " << sRef[1] << ", " << sRef[2] << ")\n";
    std::cout << "ΣF(Gen) = (" << sGen[0] << ", " << sGen[1] << ", " << sGen[2] << ")\n";
    std::cout << "max|ΔF| = " << maxDelta << "\n";
    std:: cout<<  "mean L2 = " << l2 << "\n";
    return 0;
}
