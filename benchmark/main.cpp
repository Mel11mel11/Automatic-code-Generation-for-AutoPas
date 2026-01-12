#include "Particle.h"
#include "ArrayUtils.h"
#include "ArrayMath.h"
#include "Timer.h"
#include "Functors/LJFunctorReference.h"
#include "Functors/LJFunctor_Gen.h"
#include "Functors/Soa/LJFunctor_Ref_SoA.h"
#include "Functors/Soa/LJFunctor_Gen_SoA.h"
#include "Functors/MieFunctorReference.h"
#include "Functors/MieFunctor_Gen.h"
#include "Functors/Soa/MieFunctor_Ref_SoA.h"
#include "Functors/GravFunctorReference.h"
#include "Functors/GravityFunctor_Gen.h"
#include "Functors/Soa/GravityFunctor_Ref_SoA.h"
#include "Functors/KryptonFunctorReference.h"
#include "Functors/KryptonFunctorGenerated_Gen.h"
#include "Functors/Soa/GravityFunctor_Gen_SoA.h"
#include "Functors/Soa/MieFunctor_Gen_SoA.h"
#include "Functors/Soa/KryptonFunctorGenerated_Gen_SoA.h"
#include "Functors/KryptonFunctorGenerated_Gen2.h"
#include "Functors/Soa/KryptonFunctor_Ref_SoA.h"
//------------------------------------------------------
#include "Functors/generatednew /generated_LJ_O000.hpp"
#include "Functors/generatednew /generated_LJ_O001.hpp"
#include "Functors/generatednew /generated_LJ_O010.hpp"
#include "Functors/generatednew /generated_LJ_O100.hpp"
#include "Functors/generatednew /generated_LJ_O101.hpp"
#include "Functors/generatednew /generated_LJ_O011.hpp"
#include "Functors/generatednew /generated_LJ_O110.hpp"
#include "Functors/generatednew /generated_LJ_O111.hpp"
//------------------------------------------------------
#include "Functors/generatednew /generated_Gravity_O000.hpp"
#include "Functors/generatednew /generated_Gravity_O001.hpp"
#include "Functors/generatednew /generated_Gravity_O010.hpp"
#include "Functors/generatednew /generated_Gravity_O100.hpp"
#include "Functors/generatednew /generated_Gravity_O101.hpp"
#include "Functors/generatednew /generated_Gravity_O110.hpp"
#include "Functors/generatednew /generated_Gravity_O011.hpp"
#include "Functors/generatednew /generated_Gravity_O111.hpp"
//----------------------------------------------------------

#include "Functors/generatednew /generated_Mie_O000.hpp"
#include "Functors/generatednew /generated_Mie_O001.hpp"
#include "Functors/generatednew /generated_Mie_O010.hpp"
#include "Functors/generatednew /generated_Mie_O100.hpp"
#include "Functors/generatednew /generated_Mie_O101.hpp"
#include "Functors/generatednew /generated_Mie_O110.hpp"
#include "Functors/generatednew /generated_Mie_O011.hpp"
#include "Functors/generatednew /generated_Mie_O111.hpp"
//----------------------------------------------------------
#include "Functors/generatednew /generated_Krypton_O000.hpp"
#include "Functors/generatednew /generated_Krypton_O001.hpp"
#include "Functors/generatednew /generated_Krypton_O010.hpp"
#include "Functors/generatednew /generated_Krypton_O100.hpp"
#include "Functors/generatednew /generated_Krypton_O101.hpp"
#include "Functors/generatednew /generated_Krypton_O110.hpp"
#include "Functors/generatednew /generated_Krypton_O011.hpp"
#include "Functors/generatednew /generated_Krypton_O111.hpp"
//----------------------------------------------------------




#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>
#include <string>


struct ParticleSoA {
    std::vector<double> x, y, z;
    std::vector<double> fx, fy, fz;
    std::vector<double> mass;         

    std::size_t size() const { return x.size(); }

    void resize(size_t N) {
        x.resize(N); y.resize(N); z.resize(N);
        fx.resize(N); fy.resize(N); fz.resize(N);
        mass.resize(N);             
    }
};

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

static double checksumFx(const std::vector<ParticleType>& ps)
{   double s=0; 
    for(auto& p:ps) 
    s+=std::abs(p.getF()[0]);
    return s; }

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
static void loadAoSIntoSoA(const std::vector<ParticleType> &aos, ParticleSoA &soa) {
    const size_t N = aos.size();
    soa.resize(N);

    for (size_t i = 0; i < N; ++i) {
        auto r = aos[i].getR();

        soa.x[i] = r[0];
        soa.y[i] = r[1];
        soa.z[i] = r[2];

        soa.fx[i] = 0.0;
        soa.fy[i] = 0.0;
        soa.fz[i] = 0.0;

        soa.mass[i] = aos[i].getMass(); 
    }
}


static void extractSoAIntoAoS(const ParticleSoA &soa, std::vector<ParticleType> &aos) {
    const size_t N = aos.size();
    for (size_t i = 0; i < N; ++i) {
        aos[i].setF({soa.fx[i], soa.fy[i], soa.fz[i]});
    }
}

template<class F>
static long runPairsSoA(std::vector<ParticleType>& ps, F& functor){
    ParticleSoA soa;
    soa.resize(ps.size());
    loadAoSIntoSoA(ps, soa);

    Timer t; t.start();
    functor.SoAFunctor(soa);   
    t.stop();

    extractSoAIntoAoS(soa, ps);
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

[[maybe_unused]] auto benchSoA = [&](const std::string& name, auto& functor){
    zero_all(ps);
    long t_ns = runPairsSoA(ps, functor);
    double sum = checksumFx(ps);
    std::cout << name << " (SoA)  N="<<N<<"  time="<<(t_ns/1e6)<<" ms  sumFx="<<sum << "\n";
};

[[maybe_unused]] auto bench = [&](const std::string& name, auto& functor){
    zero_all(ps); 
        zero_all(ps);

    long t_ns = runPairs(ps, functor);
    auto sums = sumAllTwoWays(ps); 
    double sum = checksumFx(ps);
    std::cout << name << "  N="<<N<<"  time="<<(t_ns/1e6)<<" ms  sumFx="<<sum << "\n";
    auto benchSoA = [&](const std::string& name, auto& functor){
    zero_all(ps);
    long t_ns = runPairsSoA(ps, functor);
    double sum = checksumFx(ps);
    std::cout << name << " (SoA)  N="<<N<<"  time="<<(t_ns/1e6)<<" ms  sumFx="<<sum << "\n";
};
 
    };

    const double cutoff = 2.5;
    const double sigma=1.0, epsilon=1.0;
    const int n=7, m=6;
    const double G=6.67430e-11;

    if (mode=="lj" || mode=="all"){
        LJFunctorReference<ParticleType> ref(sigma,epsilon, false,cutoff);
        LJFunctor_Gen<ParticleType> otolj(sigma, epsilon,false,cutoff);
        LJFunctor_Gen_SoA<ParticleSoA> ljSoa(sigma, epsilon, false,cutoff);
        LJFunctor_Ref_SoA<ParticleSoA> ljRefSoa(sigma, epsilon, false,cutoff);
        LJFunctor_Gen_O000<ParticleType> oto000(sigma, epsilon,false,cutoff);
        LJFunctor_Gen_O001<ParticleType> oto001(sigma, epsilon,false,cutoff);
        LJFunctor_Gen_O010<ParticleType> oto010(sigma, epsilon,false,cutoff);
        LJFunctor_Gen_O100<ParticleType> oto100(sigma, epsilon,false,cutoff);
        LJFunctor_Gen_O101<ParticleType> oto101(sigma, epsilon,false,cutoff);
        LJFunctor_Gen_O110<ParticleType> oto110(sigma, epsilon,false,cutoff);
        LJFunctor_Gen_O011<ParticleType> oto011(sigma, epsilon,false,cutoff);
        LJFunctor_Gen_O111<ParticleType> oto111(sigma, epsilon,false,cutoff);


        bench("LJ-REF ", ref);
     
        bench("LJ-OTO ", otolj);
        bench("LJ- no opt ", oto000);
        bench("LJ- O001 ", oto001);
        bench("LJ- O010 ", oto010);
        bench("LJ- O100 ", oto100);
        bench("LJ- O101 ", oto101);
        bench("LJ- O110 ", oto110);
        bench("LJ- O011 ", oto011);
        bench("LJ- O111 ", oto111);

        benchSoA("LJ-GEN-SOA ", ljSoa);
        benchSoA("LJ-REF-SOA ", ljRefSoa);
    }
    printf("---------------------------------------\n");
    if (mode == "mie" || mode == "all") {
        MieFunctorReference<ParticleType> mieRef(sigma, epsilon, n, m, false, cutoff);
        MieFunctor_Gen_O000<ParticleType> mieO000(sigma, epsilon, n, m, false, cutoff);
        MieFunctor_Gen_O001<ParticleType> mieO001(sigma, epsilon, n, m, false, cutoff);
        MieFunctor_Gen_O010<ParticleType> mieO010(sigma, epsilon, n, m, false, cutoff);
        MieFunctor_Gen_O100<ParticleType> mieO100(sigma, epsilon, n, m, false, cutoff);
        MieFunctor_Gen_O101<ParticleType> mieO101(sigma, epsilon, n, m, false, cutoff);
        MieFunctor_Gen_O110<ParticleType> mieO110(sigma, epsilon, n, m, false, cutoff);
        MieFunctor_Gen_O011<ParticleType> mieO011(sigma, epsilon, n, m, false, cutoff);
        MieFunctor_Gen_O111<ParticleType> mieO111(sigma, epsilon, n, m, false, cutoff);

        MieFunctor_Gen_SoA<ParticleSoA> mieSoa(sigma, epsilon, n, m, false, cutoff);
        MieFunctor_Ref_SoA<ParticleSoA> mieRefSoa(sigma, epsilon, n, m, false, cutoff);


        bench("MIE-REF ", mieRef);
        bench("MIE- no opt ", mieO000);
        bench("MIE- O001 ", mieO001);
        bench("MIE- O010 ", mieO010);
        bench("MIE- O100 ", mieO100);
        bench("MIE- O101 ", mieO101);
        bench("MIE- O110 ", mieO110);
        bench("MIE- O011 ", mieO011);
        bench("MIE- O111 ", mieO111);
        benchSoA("MIE-SOA ", mieSoa);
        benchSoA("MIE-REF-SOA ", mieRefSoa);

    }
    printf("---------------------------------------\n");
    if(mode=="krypton"|| mode=="all"){
        
    KryptonFunctorReference<ParticleType> kry_ref(1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0,
        false, cutoff);
    KryptonFunctorGenerated_Gen_O000<ParticleType> kry_000(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O001<ParticleType> kry_001(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O010<ParticleType> kry_010(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false, cutoff);

    KryptonFunctorGenerated_Gen_O100<ParticleType> kry_100(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O101<ParticleType> kry_101(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O110<ParticleType> kry_110(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O011<ParticleType> kry_011(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O111<ParticleType> kry_111(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false, cutoff);
            
    KryptonFunctorGenerated_Gen_SoA<ParticleSoA> krp_soa(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false, cutoff);
    
    KryptonFunctor_Ref_SoA<ParticleSoA> kry_ref_soa(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false, cutoff);

     bench("KRY-REF", kry_ref);
     bench("KRY- no opt", kry_000);
     bench("KRY- O001", kry_001);
     bench("KRY- O010", kry_010);
     bench("KRY- O100", kry_100);
     bench("KRY- O101", kry_101);
     bench("KRY- O110", kry_110);
     bench("KRY- O011", kry_011);
     bench("KRY- O111", kry_111);
     benchSoA("KRY-SOA", krp_soa);
     benchSoA("KRY-REF-SOA", kry_ref_soa);

    }
    printf("---------------------------------------\n");
    if (mode=="grav" || mode=="all"){
       
        GravFunctorReference<ParticleType> gRef(G,true, cutoff);
        GravityFunctor_Gen_O000<ParticleType> go000(G,true, cutoff);
        GravityFunctor_Gen_O001<ParticleType> go001(G,true, cutoff);
        GravityFunctor_Gen_O010<ParticleType> go010(G,true, cutoff);
        GravityFunctor_Gen_O100<ParticleType> go100(G,true, cutoff);
        GravityFunctor_Gen_O101<ParticleType> go101(G,true, cutoff);
        GravityFunctor_Gen_O110<ParticleType> go110(G,true, cutoff);
        GravityFunctor_Gen_O011<ParticleType> go011(G,true, cutoff);
        GravityFunctor_Gen_O111<ParticleType> go111(G,true, cutoff);
        GravityFunctor_Gen_SoA<ParticleSoA> gravSoa(G,true, cutoff);
        GravityFunctor_Ref_SoA<ParticleSoA> gravRefSoa(G,true, cutoff);

    
        bench("GRAV-REF", gRef);
        bench("GRAV- no opt", go000);
        bench("GRAV- O001", go001);
        bench("GRAV- O010", go010);
        bench("GRAV- O100", go100);
        bench("GRAV- O101", go101);
        bench("GRAV- O110", go110);
        bench("GRAV- O011", go011);
        bench("GRAV- O111", go111);
        benchSoA("GRAV-SOA", gravSoa);
        benchSoA("GRAV-REF-SOA", gravRefSoa);
    }
    // O100 simplify
    // O010 CSE
    // O001 fast_pow
    // O101 simplify + fast_pow
    // O110 simplify + CSE
    // O011 CSE + fast_pow
    // O111 full
    return 0;
}
