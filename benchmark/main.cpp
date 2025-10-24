#include "Particle.h"
#include "ArrayUtils.h"
#include "ArrayMath.h"
#include "Timer.h"

// Reference LJ
#include "Functors/LJFunctorReference.h"
#include "Functors/LJFunctorGenerated.h"
#include "Functors/MieFunctorGenerated.h"
#include "Functors/MieFunctorReference.h"
#include "Functors/GravFunctorGenerated.h"
#include "Functors/GravFunctorReference.h"

#include <vector>
#include <random>
#include <iostream>
#include <iomanip>
#include <string>

using ParticleType = Particle;

static void zero_all(std::vector<ParticleType>& ps){ for(auto& p:ps) p.setF({0,0,0}); }
static double checksumFx(const std::vector<ParticleType>& ps){ double s=0; for(auto& p:ps) s+=p.getF()[0]; return s; }

template<class F>
static long runPairs(std::vector<ParticleType>& ps, F& functor){
    Timer t; t.start();
    const size_t N=ps.size();
    for(size_t i=0;i<N;++i) for(size_t j=i+1;j<N;++j) functor.AoSFunctor(ps[i], ps[j]);
    t.stop(); return t.getTotalTime(); // ns
}

int main(int argc,char** argv){
// create a command line
    const std::string mode = (argc>=2)? argv[1] : "all"; 

    const size_t N = (argc>=3)? static_cast<size_t>(std::stoull(argv[2])) : 10000;

    
    std::vector<ParticleType> ps; ps.reserve(N);
    std::mt19937_64 rng(44);
    std::uniform_real_distribution<double> U(0.0,1.0);
    for(size_t i=0;i<N;++i) ps.emplace_back(std::array<double,3>{U(rng),U(rng),U(rng)},
                                            std::array<double,3>{0,0,0},
                                            std::array<double,3>{0,0,0},
                                            static_cast<int>(i));

    std::cout<<std::fixed<<std::setprecision(6);

    auto bench = [&](const std::string& name, auto& functor){
        zero_all(ps); // make zero all particles
        long t_ns = runPairs(ps, functor);
        double sum = checksumFx(ps);
        std::cout << name << "  N="<<N<<"  time="<<(t_ns/1e6)<<" ms  sumFx="<<sum << "\n";
    };

    
    const double sigma=1.0, epsilon=1.0;
    const int n=12, m=6;
    const double G=6.67430e-11;

    if (mode=="lj" || mode=="all"){
        LJFunctorReference<ParticleType> ref(sigma,epsilon);
        LJFunctorGenerated<ParticleType> gen(sigma,epsilon,true);
        bench("LJ-REF ", ref);
        bench("LJ-GEN ", gen);
    }
    if (mode == "mie" || mode == "all") {
        MieFunctorGenerated<ParticleType> mieFast(sigma, epsilon, n, m, true);
        MieFunctorGenerated<ParticleType> mieSafe(sigma, epsilon, n, m, true);
        MieFunctorReference<ParticleType> mieRef(sigma, epsilon, n, m, true);
    
        bench("MIE-FAST", mieFast);
        bench("MIE-SAFE", mieSafe);
        bench("MIE-REF ", mieRef);
    }
    if (mode=="grav" || mode=="all"){
        GravFunctorGenerated<ParticleType> gGen(G,true);
        GravFunctorReference<ParticleType> gRef(G,true);
        bench("GRAV-GEN", gGen);
        bench("GRAV-REF", gRef);
    }
    return 0;
}

//TODO: arguments n and m  for command line