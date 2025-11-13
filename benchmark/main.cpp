#include "Particle.h"
#include "ArrayUtils.h"
#include "ArrayMath.h"
#include "Timer.h"
#include "Functors/LJFunctorReference.h"
#include "Functors/LJFunctorGenerated.h"
#include "Functors/LJ_Oto.h"
#include "Functors/MieFunctorGenerated.h"
#include "Functors/MieFunctorReference.h"
#include "Functors/GravFunctorGenerated.h"
#include "Functors/GravFunctorReference.h"
#include "Functors/KryptonFunctorGenerated.h"
#include "Functors/KryptonFunctorReference.h"
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>
#include <string>

using ParticleType = Particle;

static void zero_all(std::vector<ParticleType>& ps){
    for (auto& p : ps) p.setF(std::array<double,3>{0.0, 0.0, 0.0});
}
static std::array<double,3> sumAllTwoWays(const std::vector<Particle>& ps){
    long double s_arr = 0.0L, s_get = 0.0L, maxAbs = 0.0L;
    size_t maxIdx = (size_t)-1, nonFinite = 0;
    for(size_t i=0;i<ps.size();++i){
        const auto  f  = ps[i].getF();   // std::array<double,3>
        const double fx_arr = f[0];
        const double fx_get = ps[i].getFx(); // 

        if(!std::isfinite(fx_arr) || !std::isfinite(fx_get)) ++nonFinite;

        s_arr += fx_arr;
        s_get += fx_get;

        long double a = std::fabsl((long double)fx_arr);
        if (a > maxAbs){ maxAbs = a; maxIdx = i; }
    }
    return { (double)s_arr, (double)s_get, (double)maxAbs };
}

static double checksumFx(const std::vector<ParticleType>& ps){ double s=0; for(auto& p:ps) s+=p.getF()[0]; return s; }

template<class F>
 
static long runPairs(std::vector<ParticleType>& ps, F& functor){ // this func wanders around the molecule pairs and they it messes the time
    Timer t; t.start();
    const size_t N = ps.size();
    for (size_t i = 0; i < N; ++i)
        for (size_t j = i + 1; j < N; ++j)
            functor.AoSFunctor(ps[i], ps[j]);
    t.stop();
    return t.getTotalTime();
}

int main(int argc,char** argv){

    const std::string mode = (argc>=2)? argv[1] : "all"; 

    const size_t N = (argc>=3)? static_cast<size_t>(std::stoull(argv[2])) : 1000;

    
    std::vector<ParticleType> ps; ps.reserve(N);
    std::mt19937_64 rng(44);
    std::uniform_real_distribution<double> U(0.0,1.0);
    for(size_t i=0;i<N;++i) ps.emplace_back(
        std::array<double,3>{U(rng),U(rng),U(rng)},
        std::array<double,3>{0,0,0},
        std::array<double,3>{0,0,0},
    static_cast<int>(i),
    1.0   
);
//
auto grav_sanity = [](){
    Particle a({0.0,0.0,0.0},{0,0,0},{0,0,0}, 0, 1.0);
    Particle b({1.0,0.0,0.0},{0,0,0},{0,0,0}, 1, 1.0);

    GravFunctorGenerated<Particle> ggen(6.67430e-11, true);  // N3 = true
    a.setF({0,0,0}); b.setF({0,0,0});
    ggen.AoSFunctor(a,b);

    auto FA = a.getF(), FB = b.getF();


    auto oldf = std::cout.flags(); auto oldp = std::cout.precision();
    std::cout.setf(std::ios::scientific); std::cout.precision(12);

    std::cout << "[SANITY] N3=true  FA=("<<FA[0]<<","<<FA[1]<<","<<FA[2]<<")  "
              << "FB=("<<FB[0]<<","<<FB[1]<<","<<FB[2]<<")  "
              << "FA+FB=("<<(FA[0]+FB[0])<<","<<(FA[1]+FB[1])<<","<<(FA[2]+FB[2])<<")\n";

    std::cout.flags(oldf); std::cout.precision(oldp);
};

// grav_sanity();  

    auto bench = [&](const std::string& name, auto& functor){
    zero_all(ps); // make zero all particles
        zero_all(ps);

    long t_ns = runPairs(ps, functor);
    auto sums = sumAllTwoWays(ps); 
    double sum = checksumFx(ps);
    std::cout << name << "  N="<<N<<"  time="<<(t_ns/1e6)<<" ms  sumFx="<<sum << "\n";
    };

    // auto sanity_one_pair = [&](){
    // ParticleType a({0.1,0.2,0.3},{0,0,0},{0,0,0}, 0);
    // ParticleType b({0.9,0.8,0.7},{0,0,0},{0,0,0}, 1);

    // LJFunctorReference<ParticleType> ref(1.0, 1.0, /*newton3*/ true);
   //  LJFunctorGenerated<ParticleType> gen(1.0, 1.0, /*newton3*/ true);

    // a.setF(std::array<double,3>{0,0,0});
    // b.setF(std::array<double,3>{0,0,0});
    // ref.AoSFunctor(a,b);
    //std::cerr<<"REF F1=("<<a.getF()[0]<<","<<a.getF()[1]<<","<<a.getF()[2]<<")  "
            // <<"F2=("<<b.getF()[0]<<","<<b.getF()[1]<<","<<b.getF()[2]<<")\n";

    // a.setF(std::array<double,3>{0,0,0});
   // b.setF(std::array<double,3>{0,0,0});
    //gen.AoSFunctor(a,b);
    //std::cerr<<"GEN F1=("<<a.getF()[0]<<","<<a.getF()[1]<<","<<a.getF()[2]<<")  "
            // <<"F2=("<<b.getF()[0]<<","<<b.getF()[1]<<","<<b.getF()[2]<<")\n";
    //}; 


    
    const double sigma=1.0, epsilon=1.0;
    const int n=7, m=6;
    const double G=6.67430e-11;
    //sanity_one_pair();
    if (mode=="lj" || mode=="all"){
        LJFunctorReference<ParticleType> ref(sigma,epsilon, true);
        LJFunctorGenerated<ParticleType> gen(sigma,epsilon,true);
        LJ_Oto<ParticleType> hadioto_lj(sigma, epsilon, true);

        bench("LJ-REF ", ref);
        bench("LJ-GEN ", gen);
        bench("LJ-OTO ", hadioto_lj);
        
    }
    if (mode == "mie" || mode == "all") {
        MieFunctorGenerated<ParticleType> mieSafe(sigma, epsilon, n, m, false);
        MieFunctorReference<ParticleType> mieRef(sigma, epsilon, n, m, false);

        bench("MIE-SAFE", mieSafe);
        bench("MIE-REF ", mieRef);
    }

    if(mode=="krypton"|| mode=="all"){
        
    KryptonFunctorReference<ParticleType> kry_gen(1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0,
        3670.0, 12600.0,  42800.0,
        false);
    KryptonFunctorGenerated<ParticleType> kryp_ref(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0,
        3670.0, 12600.0,  42800.0,
        false);
     bench("KRY-REF", kry_gen);
     bench("KRY-GEN", kryp_ref);

    }

    if (mode=="grav" || mode=="all"){
        GravFunctorGenerated<ParticleType> gGen(G,true);
        GravFunctorReference<ParticleType> gRef(G,true);
        bench("GRAV-GEN", gGen);
        bench("GRAV-REF", gRef);
    }
    return 0;
}
