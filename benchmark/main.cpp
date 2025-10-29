#include "Particle.h"
#include "ArrayUtils.h"
#include "ArrayMath.h"
#include "Timer.h"
#include <cmath>

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

// static void zero_all(std::vector<ParticleType>& ps){ for(auto& p:ps) p.setF({0,0,0}); }
static void zero_all(std::vector<ParticleType>& ps){
    for (auto& p : ps) p.setF(std::array<double,3>{0.0, 0.0, 0.0});
}

static double checksumFx(const std::vector<Particle>& ps){
    double s = 0.0;
    for (const auto& p : ps) {
        const double fx = p.getFx();  //
        if (!std::isfinite(fx)) {
            std::cerr << "[ERR] non-finite Fx: " << fx << "\n";
        }
        s += fx;
    }
    return s;
}


//static double checksumFx(const std::vector<ParticleType>& ps){ double s=0; for(auto& p:ps) s+=p.getF()[0]; return s; }

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

    const size_t N = (argc>=3)? static_cast<size_t>(std::stoull(argv[2])) : 100;

    
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
        zero_all(ps);
     std::cout << "[DBG] after zero_all sumFx=" << checksumFx(ps) << "\n";

    long t_ns = runPairs(ps, functor);
    double sum = checksumFx(ps);
    std::cout << name << "  N="<<N<<"  time="<<(t_ns/1e6)<<" ms  sumFx="<<sum << "\n";
    };

     auto sanity_one_pair = [&](){
    ParticleType a({0.1,0.2,0.3},{0,0,0},{0,0,0}, 0);
    ParticleType b({0.9,0.8,0.7},{0,0,0},{0,0,0}, 1);

    LJFunctorReference<ParticleType> ref(1.0, 1.0, /*newton3*/ true);
    LJFunctorGenerated<ParticleType> gen(1.0, 1.0, /*newton3*/ true);

    // (opsiyonel) bu fonksiyonlar varsa bayrağı göster:
    // std::cerr << std::boolalpha
    //           << "[DBG] ref.usesNewton3=" << ref.usesNewton3()
    //           << " gen.usesNewton3=" << gen.usesNewton3() << "\n";

    a.setF(std::array<double,3>{0,0,0});
    b.setF(std::array<double,3>{0,0,0});
    ref.AoSFunctor(a,b);
    std::cerr<<"REF F1=("<<a.getF()[0]<<","<<a.getF()[1]<<","<<a.getF()[2]<<")  "
             <<"F2=("<<b.getF()[0]<<","<<b.getF()[1]<<","<<b.getF()[2]<<")\n";

    a.setF(std::array<double,3>{0,0,0});
    b.setF(std::array<double,3>{0,0,0});
    gen.AoSFunctor(a,b);
    std::cerr<<"GEN F1=("<<a.getF()[0]<<","<<a.getF()[1]<<","<<a.getF()[2]<<")  "
             <<"F2=("<<b.getF()[0]<<","<<b.getF()[1]<<","<<b.getF()[2]<<")\n";
};

    
    const double sigma=1.0, epsilon=1.0;
    const int n=7, m=6;
    const double G=6.67430e-11;
    sanity_one_pair();
    if (mode=="lj" || mode=="all"){
        LJFunctorReference<ParticleType> ref(sigma,epsilon, true);
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
    if (mode=="lj_debug" || mode=="all"){ 
    const double sigma=1.0, epsilon=1.0;

    // n3=true test
    {
        LJFunctorReference<ParticleType> ref(sigma,epsilon,true);
        LJFunctorGenerated<ParticleType> gen(sigma,epsilon,true);
        //compareFunctorsOnce("[LJ n3=T]", ref, gen);
    }
    // n3=false test
    {
        LJFunctorReference<ParticleType> ref(sigma,epsilon,false);
        LJFunctorGenerated<ParticleType> gen(sigma,epsilon,false);
        //compareFunctorsOnce("[LJ n3=F]", ref, gen);
    }
    }
    if (mode=="grav" || mode=="all"){
        GravFunctorGenerated<ParticleType> gGen(G,false);
        GravFunctorReference<ParticleType> gRef(G,false);
        bench("GRAV-GEN", gGen);
        bench("GRAV-REF", gRef);
    }
    return 0;
}

//TODO: arguments n and m  for command line