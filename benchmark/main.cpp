//========================= REFERENCE FUNCTORS ==========================
#include "Functors/LJFunctorReference.h"
#include "Functors/MieFunctorReference.h"
#include "Functors/GravFunctorReference.h"
#include "Functors/KryptonFunctorReference.h"
//-----------------------SoA Reference Functors------------------------
#include "Functors/Soa/KryptonFunctor_Ref_SoA.h"
#include "Functors/Soa/MieFunctor_Ref_SoA.h"
#include "Functors/Soa/GravityFunctor_Ref_SoA.h"
#include "Functors/Soa/LJFunctor_Ref_SoA.h"
//========================LENNARD JONES OPT variations ==========================
#include "Functors/generatednew /generated_LJ_O000.hpp"
#include "Functors/generatednew /generated_LJ_O001.hpp"
#include "Functors/generatednew /generated_LJ_O010.hpp"
#include "Functors/generatednew /generated_LJ_O100.hpp"
#include "Functors/generatednew /generated_LJ_O101.hpp"
#include "Functors/generatednew /generated_LJ_O011.hpp"
#include "Functors/generatednew /generated_LJ_O110.hpp"
#include "Functors/generatednew /generated_LJ_O111.hpp"
//-------------------------CUTOFF--------------------------------
#include "Functors/without_cutoff/LJRef_wo_cutoff.h"
#include "Functors/without_cutoff_aos/LJFunctor_Gen_Opt1110.hpp"
//-------------------------CSE Optimization variants--------------
#include "Functors/new_None/generated_LJ_Opt010.hpp"
#include "Functors/cse_3/generated_LJ_Opt010.hpp"
// ----------------------SoA----------------------------
#include "Functors/generatednew /LJFunctor_Gen_O000_SoA.h"
#include "Functors/generatednew /LJFunctor_Gen_O001_SoA.h"
#include "Functors/generatednew /LJFunctor_Gen_O010_SoA.h"
#include "Functors/generatednew /LJFunctor_Gen_O011_SoA.h"
#include "Functors/generatednew /LJFunctor_Gen_O100_SoA.h"
#include "Functors/generatednew /LJFunctor_Gen_O101_SoA.h"
#include "Functors/generatednew /LJFunctor_Gen_O110_SoA.h"
#include "Functors/generatednew /LJFunctor_Gen_O111_SoA.h"
//------------------------------------------------------
//========================GRAVITY OPT variations ==========================
#include "Functors/generatednew /generated_Gravity_O000.hpp"
#include "Functors/generatednew /generated_Gravity_O001.hpp"
#include "Functors/generatednew /generated_Gravity_O010.hpp"
#include "Functors/generatednew /generated_Gravity_O100.hpp"
#include "Functors/generatednew /generated_Gravity_O101.hpp"
#include "Functors/generatednew /generated_Gravity_O110.hpp"
#include "Functors/generatednew /generated_Gravity_O011.hpp"
#include "Functors/generatednew /generated_Gravity_O111.hpp"
//-------------------------CUTOFF--------------------------------
#include "Functors/without_cutoff/GravRef_wo_cutoff.h"
#include "Functors/without_cutoff_aos/generated_Gravity_Opt1110.hpp"
//-------------------------CSE Optimization variants--------------
#include "Functors/new_None/generated_Gravity_Opt010.hpp"
#include "Functors/cse_3/generated_Gravity_Opt010.hpp"
//---------------------SoA--------------------------------
#include "Functors/generatednew /GravityFunctor_Gen_O000_SoA.h"
#include "Functors/generatednew /GravityFunctor_Gen_O001_SoA.h"
#include "Functors/generatednew /GravityFunctor_Gen_O010_SoA.h"
#include "Functors/generatednew /GravityFunctor_Gen_O011_SoA.h"
#include "Functors/generatednew /GravityFunctor_Gen_O100_SoA.h"
#include "Functors/generatednew /GravityFunctor_Gen_O101_SoA.h"
#include "Functors/generatednew /GravityFunctor_Gen_O110_SoA.h"
#include "Functors/generatednew /GravityFunctor_Gen_O111_SoA.h"
//----------------------------------------------------------
//========================MIE OPT variations ==========================
#include "Functors/generatednew /generated_Mie_O000.hpp"
#include "Functors/generatednew /generated_Mie_O001.hpp"
#include "Functors/generatednew /generated_Mie_O010.hpp"
#include "Functors/generatednew /generated_Mie_O100.hpp"
#include "Functors/generatednew /generated_Mie_O101.hpp"
#include "Functors/generatednew /generated_Mie_O110.hpp"
#include "Functors/generatednew /generated_Mie_O011.hpp"
#include "Functors/generatednew /generated_Mie_O111.hpp"
//-------------------------CUTOFF--------------------------------
#include "Functors/without_cutoff/MieRef_wo_cutoff.h"
#include "Functors/without_cutoff_aos/MieFunctor_Gen_Opt1110.hpp"
//-------------------------CSE Optimization variants--------------
#include "Functors/new_None/generated_Mie_Opt010.hpp"
#include "Functors/cse_3/generated_Mie_Opt010.hpp"

#include "Functors/generatednew /MieFunctor_Gen_O000_SoA.h" 
#include "Functors/generatednew /MieFunctor_Gen_O001_SoA.h"
#include "Functors/generatednew /MieFunctor_Gen_O010_SoA.h"
#include "Functors/generatednew /MieFunctor_Gen_O011_SoA.h"
#include "Functors/generatednew /MieFunctor_Gen_O100_SoA.h"
#include "Functors/generatednew /MieFunctor_Gen_O101_SoA.h"
#include "Functors/generatednew /MieFunctor_Gen_O110_SoA.h"
#include "Functors/generatednew /MieFunctor_Gen_O111_SoA.h"
//========================KRYPTON OPT variations==========================
#include "Functors/generatednew /generated_Krypton_O000.hpp"
#include "Functors/generatednew /generated_Krypton_O001.hpp"
#include "Functors/generatednew /generated_Krypton_O010.hpp"
#include "Functors/generatednew /generated_Krypton_O100.hpp"
#include "Functors/generatednew /generated_Krypton_O101.hpp"
#include "Functors/generatednew /generated_Krypton_O110.hpp"
#include "Functors/generatednew /generated_Krypton_O011.hpp"
#include "Functors/generatednew /generated_Krypton_O111.hpp"
//-------------------------CUTOFF--------------------------------
#include "Functors/without_cutoff/KryptonRef_wo_cutoff.h"
#include "Functors/without_cutoff_aos/generated_Krypton_Opt1110.hpp"
//-------------------------CSE Optimization variants--------------
#include "Functors/new_None/generated_Krypton_Opt010.hpp"
#include "Functors/cse_3/generated_Krypton_Opt010.hpp"
#include "Functors/Soa/KryptonFunctor_Ref_SoA.h"
#include "Functors/Soa/KryptonFunctorGenerated_Gen_SoA.h"
#include "Functors/generatednew /KryptonFunctorGenerated_Gen_O000_SoA.h"
#include "Functors/generatednew /KryptonFunctorGenerated_Gen_O001_SoA.h"
#include "Functors/generatednew /KryptonFunctorGenerated_Gen_O010_SoA.h"
#include "Functors/generatednew /KryptonFunctorGenerated_Gen_O011_SoA.h"
#include "Functors/generatednew /KryptonFunctorGenerated_Gen_O100_SoA.h"
#include "Functors/generatednew /KryptonFunctorGenerated_Gen_O101_SoA.h"
#include "Functors/generatednew /KryptonFunctorGenerated_Gen_O110_SoA.h"
#include "Functors/generatednew /KryptonFunctorGenerated_Gen_O111_SoA.h"

//=============== main Libraries and imports ============================
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <cstddef>
#include "Particle.h"
#include "ArrayUtils.h"
#include "ArrayMath.h"
#include "Timer.h"


using ParticleType = Particle;

struct ParticleSoA {
  std::vector<double> x, y, z; // positions
  std::vector<double> fx, fy, fz; // force components
  std::vector<double> mass; // mass 

  std::size_t size() const { return x.size(); }

  void resize(std::size_t N) {
    x.resize(N);
    y.resize(N);
    z.resize(N);
    fx.resize(N);
    fy.resize(N);
    fz.resize(N);
    mass.resize(N);
  }
};

struct ForceSums3D {
  std::array<double, 3> sum_getF; // get forces
  std::array<double, 3> sum_getXYZ; // get positions
  double maxAbsComponent;
  std::size_t maxIdx;
  int maxDim;  // 0=x,1=y,2=z
  std::size_t nonFinite;
};

static void zero_all(std::vector<ParticleType> &ps) { 
  // for every particle in the particle list set force to zero
  for (auto &p : ps) p.setF(std::array<double, 3>{0.0, 0.0, 0.0});
}

// Two different ways to compare the forces and the bugs

static ForceSums3D sumAllTwoWays3D(const std::vector<ParticleType> &ps) {
  long double sF[3] = {0.0L, 0.0L, 0.0L};
  long double sG[3] = {0.0L, 0.0L, 0.0L};
  long double maxAbs = 0.0L;
  std::size_t maxIdx = (std::size_t)-1;
  int maxDim = -1;
  std::size_t nonFinite = 0;

  for (std::size_t i = 0; i < ps.size(); ++i) { //  for every particle
    const auto f = ps[i].getF();  // get force components 
    const double g[3] = {ps[i].getFx(), ps[i].getFy(), ps[i].getFz()};// in other way get force components seperately

    for (int d = 0; d < 3; ++d) {
      const double fv = f[d]; // force component from getF()
      const double gv = g[d]; // force component from getFx(), getFy(), getFz()

      if (std::isfinite(fv) == false || std::isfinite(gv) == false) {
        nonFinite = nonFinite + 1;
          }
      sF[d] = sF[d] + fv;
      sG[d] = sG[d] + gv;

      const long double a = std::fabs((long double)fv);
      if (a > maxAbs) {
        maxAbs = a; // a-value is  the max
        maxIdx = i;
        maxDim = d;
      }
    }
  }

  ForceSums3D out;
  out.sum_getF = {(double)sF[0], (double)sF[1], (double)sF[2]};
  out.sum_getXYZ = {(double)sG[0], (double)sG[1], (double)sG[2]};
  out.maxAbsComponent = (double)maxAbs;
  out.maxIdx = maxIdx;
  out.maxDim = maxDim;
  out.nonFinite = nonFinite;
  return out;
}

static std::array<double, 3> checksumAbsFxyz(const std::vector<ParticleType> &ps) {
  long double sx = 0.0L, sy = 0.0L, sz = 0.0L;
  for (const auto &p : ps) {
    const auto f = p.getF();
    sx += std::abs(f[0]);
    sy += std::abs(f[1]);
    sz += std::abs(f[2]);
  }
  return {(double)sx, (double)sy, (double)sz};
}

template <class F>
// orginal particle list and its reference so that the changes can be kept 
static long runPairs(std::vector<ParticleType> &ps, F &functor) {
  Timer t;
  t.start(); // start the timer
  const std::size_t N = ps.size();
  for(std::size_t iterations = 0; iterations < 1000; iterations++)     
     for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j) functor.AoSFunctor(ps[i], ps[j]);
  t.stop();
  return t.getTotalTime()/1000;
}

static void loadAoSIntoSoA(const std::vector<ParticleType> &aos, ParticleSoA &soa) {
  const std::size_t N = aos.size();
  soa.resize(N);
  for (std::size_t i = 0; i < N; ++i) {
    const auto r = aos[i].getR();
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
  const std::size_t N = aos.size();
  for (std::size_t i = 0; i < N; ++i) {
    aos[i].setF({soa.fx[i], soa.fy[i], soa.fz[i]});
  }
}

template <class F>
static long runPairsSoA(std::vector<ParticleType> &ps, F &functor) {
  ParticleSoA soa;
  loadAoSIntoSoA(ps, soa);

  Timer t;
  t.start();
  for (std::size_t iter = 0; iter < 1000; ++iter) {
    std::fill(soa.fx.begin(), soa.fx.end(), 0.0); // set zero
    std::fill(soa.fy.begin(), soa.fy.end(), 0.0);
    std::fill(soa.fz.begin(), soa.fz.end(), 0.0);
    functor.SoAFunctor(soa);
  }
  t.stop();

  extractSoAIntoAoS(soa, ps);
  return t.getTotalTime()/1000;
}

// =======================
// AoS vs SoA verification
// =======================
struct ForceDiffReport {
  double maxAbsDiff = 0.0;
  double rmsDiff = 0.0;
  std::size_t maxIdx = (std::size_t)-1;
  int maxDim = -1;
  std::array<double, 3> sumA{0.0, 0.0, 0.0};
  std::array<double, 3> sumB{0.0, 0.0, 0.0};
  std::size_t nonFinite = 0;
};

static ForceDiffReport compareForces(const std::vector<ParticleType> &A,
                                    const std::vector<ParticleType> &B) {
  ForceDiffReport rep;
  long double sse = 0.0L;
  long double cnt = 0.0L;

  for (std::size_t i = 0; i < A.size(); ++i) {
    const auto fa = A[i].getF();
    const auto fb = B[i].getF();

    for (int d = 0; d < 3; ++d) {
      const double a = fa[d];
      const double b = fb[d];

      if (!std::isfinite(a) || !std::isfinite(b)) ++rep.nonFinite;

      rep.sumA[d] += a;
      rep.sumB[d] += b;

      const double diff = a - b;
      const double ad = std::abs(diff);

      if (ad > rep.maxAbsDiff) {
        rep.maxAbsDiff = ad;
        rep.maxIdx = i;
        rep.maxDim = d;
      }

      sse += (long double)diff * (long double)diff;
      cnt += 1.0L;
    }
  }

  rep.rmsDiff = (cnt > 0.0L) ? std::sqrt((double)(sse / cnt)) : 0.0;
  return rep;
}

template <class FAOS, class FSOA>
static void verifyAoSvsSoA(const std::string &name,
                           const std::vector<ParticleType> &ps0,
                           FAOS &aosFunctor,
                           FSOA &soaFunctor,
                           double atol = 1e-10,
                           double rtol = 1e-10) {
  std::vector<ParticleType> psA = ps0;
  std::vector<ParticleType> psB = ps0;

  zero_all(psA);
  runPairs(psA, aosFunctor);

  zero_all(psB);
  runPairsSoA(psB, soaFunctor);

  const auto rep = compareForces(psA, psB);

  // tolerance based on reference magnitude
  double maxAbsF = 0.0;
  for (const auto &p : psA) {
    const auto f = p.getF();
    maxAbsF = std::max(maxAbsF, std::abs(f[0]));
    maxAbsF = std::max(maxAbsF, std::abs(f[1]));
    maxAbsF = std::max(maxAbsF, std::abs(f[2]));
  }
  const double tol = atol + rtol * maxAbsF;

  std::cout << "[VERIFY] " << name
            << "  maxAbsDiff=" << std::setprecision(12) << rep.maxAbsDiff
            << "  rmsDiff=" << rep.rmsDiff
            << "  tol=" << tol
            << "  sumF_A=(" << rep.sumA[0] << "," << rep.sumA[1] << "," << rep.sumA[2] << ")"
            << "  sumF_B=(" << rep.sumB[0] << "," << rep.sumB[1] << "," << rep.sumB[2] << ")";

  if (rep.nonFinite) std::cout << "  nonFinite=" << rep.nonFinite;

  if (rep.maxAbsDiff > tol) {
    std::cout << "  ==> FAIL at i=" << rep.maxIdx << " dim=" << rep.maxDim << "\n";
  } else {
    std::cout << "  ==> PASS\n";
  }
}
// start of the main
int main(int argc, char **argv) {
  const std::string mode = (argc >= 2) ? argv[1] : "all"; // for different modes lj, mie, krypton, grav
  const std::size_t N = (argc >= 3) ? static_cast<std::size_t>(std::stoull(argv[2])) : 1000; // input particle number

  std::vector<ParticleType> ps; // particle list
  ps.reserve(N); // reserve space for N particles

  std::mt19937_64 rng(44); // random 
  std::uniform_real_distribution<double> U(0.0, 10.0); // cube with [0,1]
//Particles are initialized at uniformly random positions within a
// unit cube rather than on a regular grid to avoid introducing artificial spatial ordering effects.
  for (std::size_t i = 0; i < N; ++i) {
    ps.emplace_back( // place the particle at random position
        std::array<double, 3>{U(rng), U(rng), U(rng)}, // position
        std::array<double, 3>{0.0, 0.0, 0.0}, // velocity
        std::array<double, 3>{0.0, 0.0, 0.0}, // force
        static_cast<int>(i), // particle id
        1.0); // mass
  }

  // snapshot for verification: same initial positions, clean state
  const std::vector<ParticleType> ps0 = ps; // copy of the original particle list

  [[maybe_unused]] auto benchSoA = [&](const std::string &name, auto &functor) {
    std::vector<ParticleType> tmp = ps0;
    zero_all(tmp);
    long t_ns = runPairsSoA(tmp, functor);
    const auto sumAbs = checksumAbsFxyz(tmp);
    std::cout << name << " (SoA)  N=" << N << "  time=" << (t_ns / 1e6) << " ms"
              << "  sum|F|=(" << sumAbs[0] << "," << sumAbs[1] << "," << sumAbs[2] << ")\n";
  };

  [[maybe_unused]] auto bench = [&](const std::string &name, auto &functor) {
    std::vector<ParticleType> tmp = ps0;
    zero_all(tmp);
    long t_ns = runPairs(tmp, functor);

    const auto sums3 = sumAllTwoWays3D(tmp);
    const auto sumAbs = checksumAbsFxyz(tmp);

    std::cout << name << "  N=" << N << "  time=" << (t_ns / 1e6) << " ms"
              << "  sum|F|=(" << sumAbs[0] << "," << sumAbs[1] << "," << sumAbs[2] << ")";
    if (sums3.nonFinite) std::cout << "  nonFinite=" << sums3.nonFinite;
    std::cout << "\n";
  };

  const double cutoff = 2.0;
  const double sigma = 1.0, epsilon = 1.0;
  const int n = 7, m = 6;
  const double G = 6.67430e-11;
  const bool newton3 = false;

  if (mode == "lj" || mode == "all") {
    LJFunctorReference<ParticleType> ref(sigma, epsilon, newton3, cutoff);
    LJFunctorReference_wo_cutoff<ParticleType> ref_wo(sigma, epsilon, newton3);
    LJFunctor_Gen_Opt1110_without<ParticleType> ref_wo_opt(sigma, epsilon, newton3);
    LJFunctor_Gen_O000<ParticleType> oto000(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O001<ParticleType> oto001(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O010<ParticleType> oto010(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_Opt010_None <ParticleType> oto010_none(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O100<ParticleType> oto100(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O101<ParticleType> oto101(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O110<ParticleType> oto110(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O011<ParticleType> oto011(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O111<ParticleType> oto111(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_Opt010_powsimp<ParticleType> oto010_powsimp(sigma, epsilon, newton3, cutoff);

    LJFunctor_Ref_SoA<ParticleSoA> ljRefSoa(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O000_SoA<ParticleSoA> lj_o000_soa(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O001_SoA<ParticleSoA> lj_o001_soa(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O010_SoA<ParticleSoA> lj_o010_soa(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O100_SoA<ParticleSoA> lj_o100_soa(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O101_SoA<ParticleSoA> lj_o101_soa(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O110_SoA<ParticleSoA> lj_o110_soa(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O011_SoA<ParticleSoA> lj_o011_soa(sigma, epsilon, newton3, cutoff);
    LJFunctor_Gen_O111_SoA<ParticleSoA> lj_o111_soa(sigma, epsilon, newton3, cutoff);
    
    
    bench("LJ-REF with cutoff reference",    ref);
     std::cout << "---------------------------------------\n";
    bench("LJ-REF without cutoff reference", ref_wo);
     std::cout << "---------------------------------------\n";
    bench("LJ-REF without cutoff 0111",      ref_wo_opt);
     std::cout << "---------------------------------------\n";
    bench("LJ-O111 with cutoff 0111",        oto111);
     std::cout << "---------------------------------------\n";
    bench("LJ-O000 (no opt)",                oto000);
     std::cout << "---------------------------------------\n";
    bench("LJ-O001 (fast_pow)",              oto001);
     std::cout << "---------------------------------------\n";
    bench("LJ-O010 (CSE)",                   oto010);
     std::cout << "---------------------------------------\n";
    bench("LJ-O010 (CSE+none)",              oto010_none);
     std::cout << "---------------------------------------\n";
    bench("LJ-O010 (CSE+powsimp)",           oto010_powsimp);
     std::cout << "---------------------------------------\n";
    bench("LJ-O100 (simplify)",              oto100);
     std::cout << "---------------------------------------\n";
    bench("LJ-O101 (simplify+fast_pow)",     oto101);
     std::cout << "---------------------------------------\n";
    bench("LJ-O110 (simplify+CSE)",          oto110);
     std::cout << "---------------------------------------\n";
    bench("LJ-O011 (CSE+fast_pow)",          oto011);
     std::cout << "---------------------------------------\n";

  std::cout << "🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢\n";

    benchSoA("LJ-REF-SOA (no opt)", ljRefSoa);
    benchSoA("LJ-O000-SOA (no opt)", lj_o000_soa);
    benchSoA("LJ-O001-SOA (fast_pow)", lj_o001_soa);
    benchSoA("LJ-O010-SOA (CSE)", lj_o010_soa);
    benchSoA("LJ-O100-SOA (simplify)", lj_o100_soa);
    benchSoA("LJ-O101-SOA (simplify+fast_pow)", lj_o101_soa);
    benchSoA("LJ-O110-SOA (simplify+CSE)", lj_o110_soa);
    benchSoA("LJ-O011-SOA (CSE+fast_pow)", lj_o011_soa);
    benchSoA("LJ-O111-SOA (full)", lj_o111_soa);

  }

  std::cout << "---------------------------------------\n";

  if (mode == "mie" || mode == "all") {
    MieFunctorReference<ParticleType> mieRef(sigma, epsilon, n, m, false, cutoff);
    MieFunctorReference_wo_cutoff<ParticleType> mieRef_wo(sigma, epsilon, n, m, false);
    MieFunctor_Gen_O000<ParticleType> mieO000(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O001<ParticleType> mieO001(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O010<ParticleType> mieO010(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_Opt010_None<ParticleType> mieO010_none(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O100<ParticleType> mieO100(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O101<ParticleType> mieO101(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O110<ParticleType> mieO110(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O011<ParticleType> mieO011(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O111<ParticleType> mieO111(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_Opt010_powsimp<ParticleType> mieO010_powsimp(sigma, epsilon, n, m, false, cutoff);

    MieFunctor_Ref_SoA<ParticleSoA> mieRefSoa(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_Opt1110_without<ParticleType> mie_wo_opt(sigma, epsilon, n, m, false);
    MieFunctor_Gen_O000_SoA<ParticleSoA> mie_o000_soa(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O001_SoA<ParticleSoA> mie_o001_soa(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O010_SoA<ParticleSoA> mie_o010_soa(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O100_SoA<ParticleSoA> mie_o100_soa(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O101_SoA<ParticleSoA> mie_o101_soa(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O110_SoA<ParticleSoA> mie_o110_soa(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O011_SoA<ParticleSoA> mie_o011_soa(sigma, epsilon, n, m, false, cutoff);
    MieFunctor_Gen_O111_SoA<ParticleSoA> mie_o111_soa(sigma, epsilon, n, m, false, cutoff);
    
    bench("MIE-REF with cutoff reference", mieRef);
    std::cout << "---------------------------------------\n";
    bench("MIE-REF without cutoff reference", mieRef_wo);
    std::cout << "---------------------------------------\n";
    bench("MIE-REF 0111 without cutoff ", mie_wo_opt);
    std::cout << "---------------------------------------\n";
    bench("MIE-O111 with cutoff ", mieO111);
     std::cout << "---------------------------------------\n";
    bench("MIE-O000 (no opt)", mieO000);
    std::cout << "---------------------------------------\n";
    bench("MIE-O001 (fast_pow)", mieO001);
    std::cout << "---------------------------------------\n";
    bench("MIE-O010 (CSE)", mieO010);
    std::cout << "---------------------------------------\n";
    bench("MIE-O010_None (CSE_None)", mieO010_none); 
    std::cout << "---------------------------------------\n";
    bench("MIE-O010_powsimp (CSE_powsimp)", mieO010_powsimp); 
    std::cout << "---------------------------------------\n";
    bench("MIE-O100 (simplify)", mieO100);
    std::cout << "---------------------------------------\n";
    bench("MIE-O101(simplify+fast_pow)", mieO101);
    std::cout << "---------------------------------------\n";
    bench("MIE-O110 (simplify+CSE)", mieO110);
    std::cout << "---------------------------------------\n";
    bench("MIE-O011 (CSE+fast_pow)", mieO011);
   
  std::cout << "🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢\n";

    benchSoA("MIE-REF-SOA", mieRefSoa);
    benchSoA("MIE-O000-SOA (no opt)", mie_o000_soa);
    benchSoA("MIE-O001-SOA (fast_pow)", mie_o001_soa);
    benchSoA("MIE-O010-SOA (CSE)", mie_o010_soa);
  
    benchSoA("MIE-O100-SOA (simplify)", mie_o100_soa);
    benchSoA("MIE-O101-SOA(simplify+fast_pow)", mie_o101_soa);
    benchSoA("MIE-O110-SOA (simplify+CSE)", mie_o110_soa);
    benchSoA("MIE-O011-SOA (CSE+fast_pow)", mie_o011_soa);
    benchSoA("MIE-O111-SOA (full)", mie_o111_soa);

  }

  std::cout << "---------------------------------------\n";

  if (mode == "krypton" || mode == "all") {
    KryptonFunctorReference<ParticleType> kry_ref(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                 64.3, 307.2, 1096.0, false, cutoff);

    KryptonFunctorReference_wo_cutoff<ParticleType> kry_ref_wo(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                         64.3, 307.2, 1096.0, false);
    KryptonFunctorGenerated_Gen_Opt1110_without<ParticleType> kry_gen_full_opt(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false);
    KryptonFunctorGenerated_Gen_O000<ParticleType> kry_000(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O001<ParticleType> kry_001(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O010<ParticleType> kry_010(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false, cutoff);

    KryptonFunctorGenerated_Gen_Opt010_None<ParticleType> kry_010_none(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false, cutoff);

    KryptonFunctorGenerated_Gen_Opt010_powsimp<ParticleType> kry_010_powsimp(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false, cutoff);

    KryptonFunctorGenerated_Gen_O100<ParticleType> kry_100(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O101<ParticleType> kry_101(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O110<ParticleType> kry_110(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O011<ParticleType> kry_011(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O111<ParticleType> kry_111(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                          64.3, 307.2, 1096.0, false, cutoff);

    KryptonFunctorGenerated_Gen_SoA<ParticleSoA> krp_soa(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                        64.3, 307.2, 1096.0, false, cutoff);

    KryptonFunctor_Ref_SoA<ParticleSoA> kry_ref_soa(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                   64.3, 307.2, 1096.0, false, cutoff);
    
    KryptonFunctorGenerated_Gen_O000_SoA<ParticleSoA> kry_o000_soa(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                                   64.3, 307.2, 1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O001_SoA<ParticleSoA> kry_o001_soa(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                                     64.3, 307.2, 1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O010_SoA<ParticleSoA> kry_o010_soa(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                                        64.3, 307.2, 1096.0, false, cutoff);
    KryptonFunctorGenerated_Gen_O100_SoA<ParticleSoA> kry_o100_soa(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                                        64.3, 307.2, 1096.0, false, cutoff);    
    KryptonFunctorGenerated_Gen_O101_SoA<ParticleSoA> kry_o101_soa(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                                        64.3, 307.2, 1096.0, false, cutoff);                

    KryptonFunctorGenerated_Gen_O110_SoA<ParticleSoA> kry_o110_soa(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                                        64.3, 307.2, 1096.0, false, cutoff);    
    KryptonFunctorGenerated_Gen_O011_SoA<ParticleSoA> kry_o011_soa(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                                        64.3, 307.2, 1096.0, false, cutoff);    
    KryptonFunctorGenerated_Gen_O111_SoA<ParticleSoA> kry_o111_soa(1.213e4, 2.821, -0.748, 0.972, 13.29,
                                                                        64.3, 307.2, 1096.0, false, cutoff);

    bench("KRY-REF with cutoff reference", kry_ref);
        std::cout << "---------------------------------------\n";
    bench("KRY-REF without cutoff reference", kry_ref_wo);
        std::cout << "---------------------------------------\n";
    bench("KRY- full without cutoff", kry_gen_full_opt);
        std::cout << "---------------------------------------\n";
    bench("KRY-O111 with cutoff (full)", kry_111);
       std::cout << "---------------------------------------\n";
    bench("KRY-O000 (no opt)", kry_000);
        std::cout << "---------------------------------------\n";
    bench("KRY-O001 (fast_pow)", kry_001);
        std::cout << "---------------------------------------\n";
    bench("KRY-O010 (CSE)", kry_010);
        std::cout << "---------------------------------------\n";
    bench("KRY-O010 (CSE) none", kry_010_none);
        std::cout << "---------------------------------------\n";
    bench("KRY-O010 (CSE) powsimp", kry_010_powsimp);
        std::cout << "---------------------------------------\n";
    bench("KRY-O100 (simplify)", kry_100);
        std::cout << "---------------------------------------\n";
    bench("KRY-O101(simplify+fast_pow)", kry_101);
        std::cout << "---------------------------------------\n";
    bench("KRY-O110", kry_110);
        std::cout << "---------------------------------------\n";
    bench("KRY-O011 (CSE+fast_pow)", kry_011);
        std::cout << "---------------------------------------\n";

        std::cout << "---------------------------------------\n";
  std::cout << "🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢\n";

    benchSoA("KRY-REF-SOA", kry_ref_soa);
    benchSoA("KRY-O000-SOA (no opt)", kry_o000_soa);
    benchSoA("KRY-O001-SOA (fast_pow)", kry_o001_soa);
    benchSoA("KRY-O010-SOA (CSE)", kry_o010_soa);
    benchSoA("KRY-O100-SOA (simplify)", kry_o100_soa);
    benchSoA("KRY-O101-SOA (simplify+fast_pow)", kry_o101_soa);
    benchSoA("KRY-O110-SOA (simplify+CSE)", kry_o110_soa);
    benchSoA("KRY-O011-SOA (CSE+fast_pow)", kry_o011_soa);
    benchSoA("KRY-O111-SOA (full)", kry_o111_soa);
   
  }

  std::cout << "---------------------------------------\n";

  if (mode == "grav" || mode == "all") {
    GravFunctorReference<ParticleType> gRef(G, true, cutoff);
    GravFunctorReference_wo_cutoff<ParticleType> gRef_wo(G, true);
    GravityFunctor_Gen_Opt1110_without<ParticleType> gRef_wo_opt(G, true);
    GravityFunctor_Gen_O000<ParticleType> go000(G, true, cutoff);
    GravityFunctor_Gen_O001<ParticleType> go001(G, true, cutoff);
    GravityFunctor_Gen_Opt010_None<ParticleType> go001_none(G, true, cutoff);
    GravityFunctor_Gen_Opt010_powsimp<ParticleType> go001_powsimp(G, true, cutoff);
    GravityFunctor_Gen_O010<ParticleType> go010(G, true, cutoff);
    GravityFunctor_Gen_O100<ParticleType> go100(G, true, cutoff);
    GravityFunctor_Gen_O101<ParticleType> go101(G, true, cutoff);
    GravityFunctor_Gen_O110<ParticleType> go110(G, true, cutoff);
    GravityFunctor_Gen_O011<ParticleType> go011(G, true, cutoff);
    GravityFunctor_Gen_O111<ParticleType> go111(G, true, cutoff);
    

    GravityFunctor_Ref_SoA<ParticleSoA> gravRefSoa(G, true, cutoff);
    GravityFunctor_Gen_O000_SoA<ParticleSoA> grav_o000_soa(G, true, cutoff);
    GravityFunctor_Gen_O001_SoA<ParticleSoA> grav_o001_soa(G, true, cutoff);
    GravityFunctor_Gen_O010_SoA<ParticleSoA> grav_o010_soa(G, true, cutoff);
    GravityFunctor_Gen_O100_SoA<ParticleSoA> grav_o100_soa(G, true, cutoff);
    GravityFunctor_Gen_O101_SoA<ParticleSoA> grav_o101_soa(G, true, cutoff);
    GravityFunctor_Gen_O110_SoA<ParticleSoA> grav_o110_soa(G, true, cutoff);
    GravityFunctor_Gen_O011_SoA<ParticleSoA> grav_o011_soa(G, true, cutoff);
    GravityFunctor_Gen_O111_SoA<ParticleSoA> grav_o111_soa(G, true, cutoff);


  
    // Bench AoS
    bench("GRAV-REF  with cutoff reference", gRef);
    std::cout << "---------------------------------------\n";
    bench("GRAV-REF  without cutoff reference", gRef_wo);
        std::cout << "---------------------------------------\n";
    bench("GRAV-0111 without cutoff generated", gRef_wo_opt);
        std::cout << "---------------------------------------\n";
    bench("GRAV-O111  with cutoff generated", go111);
    std::cout << "---------------------------------------\n";
    bench("GRAV-O000 (no opt)", go000);
        std::cout << "---------------------------------------\n";
    bench("GRAV-O001 (fast_pow)", go001);
        std::cout << "---------------------------------------\n";
    bench("GRAV-O010 (CSE)", go010);
        std::cout << "---------------------------------------\n";
    bench("GRAV-O010 (CSE) none", go001_none);
        std::cout << "---------------------------------------\n";
    bench("GRAV-O010 (CSE) powsimp", go001_powsimp);
        std::cout << "---------------------------------------\n";
    bench("GRAV-O100 (simplify)", go100);
        std::cout << "---------------------------------------\n";
    bench("GRAV-O101 (simplify+fast_pow)", go101);
        std::cout << "---------------------------------------\n";
    bench("GRAV-O110 (simplify+CSE)", go110);
        std::cout << "---------------------------------------\n";
    bench("GRAV-O011 (CSE+fast_pow)", go011);
  std::cout << "🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢🟢\n";
    // Bench SoA
   // benchSoA("GRAV-GEN-SOA", gravSoa);
    benchSoA("GRAV-O000-SOA (no opt)", grav_o000_soa);
    benchSoA("GRAV-O001-SOA (fast_pow)", grav_o001_soa);
    benchSoA("GRAV-O010-SOA (CSE)", grav_o010_soa);
    benchSoA("GRAV-O100-SOA (simplify)", grav_o100_soa);
    benchSoA("GRAV-O101-SOA (simplify+fast_pow)", grav_o101_soa);
    benchSoA("GRAV-O110-SOA (simplify+CSE)", grav_o110_soa);
    benchSoA("GRAV-O011-SOA (CSE+fast_pow)", grav_o011_soa);
    benchSoA("GRAV-O111-SOA (full)", grav_o111_soa);

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
