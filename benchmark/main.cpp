#include "Particle.h"
#include "ArrayUtils.h"
#include "ArrayMath.h"
#include "Timer.h"
#include "Functors/LJFunctorReference.h"
#include "Functors/LJFunctorGenerated.h"
#include "Functors/LJFunctor_Gen.h"
#include "Functors/Soa/LJFunctor_Ref_SoA.h"
#include "Functors/Soa/LJFunctor_Gen_SoA.h"
#include "Functors/MieFunctorGenerated.h"
#include "Functors/MieFunctorReference.h"
#include "Functors/MieFunctor_Gen.h"
#include "Functors/Soa/MieFunctor_Ref_SoA.h"
#include "Functors/GravFunctorGenerated.h"
#include "Functors/GravFunctorReference.h"
#include "Functors/GravityFunctor_Gen.h"
#include "Functors/Soa/GravityFunctor_Ref_SoA.h"
#include "Functors/KryptonFunctorGenerated.h"
#include "Functors/KryptonFunctorReference.h"
#include "Functors/KryptonFunctorGenerated_Gen.h"
#include "Functors/Soa/GravityFunctor_Gen_SoA.h"
#include "Functors/Soa/MieFunctor_Gen_SoA.h"
#include "Functors/Soa/KryptonFunctorGenerated_Gen_SoA.h"
#include "Functors/KryptonFunctorGenerated_Gen2.h"
#include "Functors/Soa/KryptonFunctor_Ref_SoA.h"
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
    
    const double sigma=1.0, epsilon=1.0;
    const int n=7, m=6;
    const double G=6.67430e-11;

    if (mode=="lj" || mode=="all"){
        LJFunctorReference<ParticleType> ref(sigma,epsilon, false);
        LJFunctorGenerated<ParticleType> gen(sigma,epsilon, false);
        LJFunctor_Gen<ParticleType> otolj(sigma, epsilon,false);
        LJFunctor_Gen_SoA<ParticleSoA> ljSoa(sigma, epsilon, false);
        LJFunctor_Ref_SoA<ParticleSoA> ljRefSoa(sigma, epsilon, false);

        bench("LJ-REF ", ref);
        bench("LJ-GEN ", gen);
        bench("LJ-OTO ", otolj);
        benchSoA("LJ-GEN-SOA ", ljSoa);
        benchSoA("LJ-REF-SOA ", ljRefSoa);
    }

    if (mode == "mie" || mode == "all") {
        MieFunctorGenerated<ParticleType> mieSafe(sigma, epsilon, n, m, false);
        MieFunctorReference<ParticleType> mieRef(sigma, epsilon, n, m, false);
        MieFunctor_Gen<ParticleType> mieOto(sigma, epsilon, n, m, false);
        MieFunctor_Gen_SoA<ParticleSoA> mieSoa(sigma, epsilon, n, m, false);
        MieFunctor_Ref_SoA<ParticleSoA> mieRefSoa(sigma, epsilon, n, m, false);

        bench("MIE-SAFE", mieSafe);
        bench("MIE-REF ", mieRef);
        bench("MIE-Oto ", mieOto);
        benchSoA("MIE-SOA ", mieSoa);
        benchSoA("MIE-REF-SOA ", mieRefSoa);

    }

    if(mode=="krypton"|| mode=="all"){
        
    KryptonFunctorReference<ParticleType> kry_gen(1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0,
        false);
    KryptonFunctorGenerated<ParticleType> kryp_ref(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false);

    KryptonFunctorGenerated_Gen<ParticleType> krp_oto(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0,
        false);  // sadece bool 
    KryptonFunctorGenerated_Gen_SoA<ParticleSoA> krp_soa(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false);
    KryptonFunctorGenerated_Gen2<ParticleType> kry_gen2(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0,
        false);  //    
    KryptonFunctor_Ref_SoA<ParticleSoA> kry_ref_soa(
        1.213e4, 2.821, -0.748, 0.972, 13.29,
        64.3,  307.2,  1096.0, false);

     bench("KRY-REF", kry_gen);
     bench("KRY-GEN", kryp_ref);
     bench("KRY-AUTO", krp_oto);
     bench("KRY-GEN2", kry_gen2);
     benchSoA("KRY-SOA", krp_soa);
     benchSoA("KRY-REF-SOA", kry_ref_soa);

    }

    if (mode=="grav" || mode=="all"){
        GravFunctorGenerated<ParticleType> gGen(G,true);
        GravFunctorReference<ParticleType> gRef(G,true);
        GravityFunctor_Gen<ParticleType> go(G,true);
        GravityFunctor_Gen_SoA<ParticleSoA> gravSoa(G,true);
        GravityFunctor_Ref_SoA<ParticleSoA> gravRefSoa(G,true);

        bench("GRAV-GEN", gGen);
        bench("GRAV-REF", gRef);
        bench("GRAV-OTO", go);
        benchSoA("GRAV-SOA", gravSoa);
        benchSoA("GRAV-REF-SOA", gravRefSoa);
    }
    return 0;
}
