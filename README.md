# Automatic Code Generation for AutoPas

A Python-based code generator that takes potential energy function definitions from a YAML file, symbolically derives the force expression, applies configurable optimizations, and writes ready-to-compile C++ functor headers for use inside the [AutoPas](https://github.com/AutoPas/AutoPas) particle simulation library.

---

## Project Structure

```
.
├── Generalization/
│   ├── main.py            # Entry point: reads YAML, runs codegen for all potentials × all opts
│   ├── yaml_loader.py     # Loads single or multi-document YAML potential specs
│   ├── validation.py      # Validates and normalizes each potential definition
│   ├── emit_header.py     # Assembles the final C++ AoS functor header as a string
│   ├── replace.py         # Post-processing: adds std:: namespace, replaces pow with fast_pow
│   ├── test.yaml          # Example YAML with LJ, Gravity, Mie, and Krypton potentials
│   └── templates/
│       └── soa_functor.h.mako   # Mako template for Structure-of-Arrays (SoA) functors
│
└── benchmark/
    ├── main.cpp           # Benchmark runner: includes all generated headers and times them
    ├── Particle.h         # Minimal particle class (position, force, mass)
    ├── Timer.h            # High-resolution wall-clock timer
    ├── ArrayMath.h        # Inline vector math helpers
    ├── Makefile           # Build script for the benchmark binary
    └── Functors/
        ├── Functor.h               # Abstract base class for all functors
        ├── FastPow.hpp             # fast_pow: repeated-squaring for small integer exponents
        ├── *FunctorReference.h     # Hand-written reference implementations
        ├── Soa/                    # Hand-written SoA reference functors
        ├── generatednew/           # Generated AoS + SoA headers (output of codegen)
        ├── without_cutoff/         # Reference functors with cutoff guard removed
        └── without_cutoff_aos/     # Generated AoS functors without cutoff (for comparison)
```

---

## How It Works

1. **YAML spec** — each potential is described by a symbolic energy expression `U(r)` and a list of parameters.
2. **Symbolic derivation** — `main.py` uses [SymPy](https://www.sympy.org/) to compute the force `F = -dU/dr` analytically.
3. **Optimization pass** — one or more of three transformations are applied:
   - `simplify_elim` — simplify the expression after substituting `r → 1/inv_r`
   - `cse` — common subexpression elimination (introduces `const double x0 = ...` temporaries)
   - `fast_pow` — replace `std::pow(base, n)` with `fast_pow(base, n)` for small integer `n`
4. **Codegen** — `emit_header.py` assembles a complete C++ class (AoS functor). If `soa: true` is set in the YAML, a second header is generated from the Mako template.
5. **Benchmark** — `benchmark/main.cpp` includes all generated headers and measures wall-clock time per functor across a random particle pair list.

All eight optimization combinations (`Opt000` through `Opt111`) are generated in a single run, making it easy to compare their performance in the benchmark.

---

## Supported Potentials

| Name | Expression | Special features |
|------|-----------|-----------------|
| Lennard-Jones | `4ε((σ/r)^12 − (σ/r)^6)` | — |
| Mie | `C·ε((σ/r)^n − (σ/r)^m)` | prefactor `C` computed at construction |
| Gravity | `−G·m₁·m₂ / r` | reads particle masses from `getMass()` |
| Krypton (Tang–Toennies) | short-range repulsion + damped dispersion series | truncated Taylor series `TT(n, b·r)`, C12/C14/C16 derived from C6/C8/C10 |

---

## Usage

### 1. Generate headers

```bash
cd Generalization
python3 main.py test.yaml ../benchmark/Functors/generatednew/
```

This creates one `.hpp`/`.h` file per potential per optimization variant in the output directory.

### 2. Build and run the benchmark

```bash
cd benchmark
make
./benchmark
```

---

## YAML Spec Format

```yaml
name: LennardJones
expression: "4*epsilon*((sigma/r)**12 - (sigma/r)**6)"
parameters:
  sigma: 1.0
  epsilon: 1.0
output:
  classname: LJFunctor_Gen
  filename: generated_LJ.hpp
options:
  soa: false           # set to true to also generate a SoA header
  newton3: true        # default Newton's third law flag in constructor
  cutoff: 2.5          # optional cutoff radius (omit to disable)
  avoid_r2_zero: 1e-24 # minimum r^2 guard against division by zero
```

Multiple potentials can be placed in a single file separated by `---`.

---

## Dependencies

- Python 3.10+
- [SymPy](https://www.sympy.org/) — symbolic math
- [PyYAML](https://pyyaml.org/) — YAML parsing
- [Mako](https://www.makotemplates.org/) — templating for SoA headers
- A C++17 compiler (GCC or Clang) for building the benchmark
