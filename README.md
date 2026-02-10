<!-- need to use divs to center images -->
<div align="center">
  <img src="./docs/images/planck.png">
</div>

### Hartree-Fock

<p align="justify"> Hartree-Fock is a from-scratch C++ implementation of the Hartree–Fock (HF) method for ab initio electronic structure calculations. The project is primarily designed as an educational and research-oriented codebase, focusing on clarity, modern C++ design, and performance-conscious implementation rather than feature parity with large production quantum chemistry packages.</p>

<p align="justify"> This repository represents a complete rewrite of the earlier project <a href="https://github.com/HemanthHaridas/plank.py">plank.py</a>, which combined CPython and Cython. While functional, the earlier version suffered from significant performance bottlenecks and architectural limitations. The current version addresses these issues through: </p>

* A pure C++ implementation
* Modern C++ (C++23) design choices
* OpenMP parallelism in computationally intensive regions
* Cleaner separation of concerns (I/O, basis sets, integrals, symmetry, etc.)

<p align="justify" style="font-weight: bold; font-style: italic"> Hartree-Fock is not intended to replace established quantum chemistry packages such as Gaussian, GAMESS, or NWChem. Its goals are transparency, hackability, and pedagogical value. </p>

#### Features
* Restricted Hartree–Fock (RHF)
* Gaussian-type orbital (GTO) basis sets (e.g., STO-3G)
* Analytical one- and two-electron integrals
* OpenMP parallelization for improved scalability
* Simple, human-readable input format

Planned and experimental features include:

* Unrestricted Hartree–Fock (UHF)
* Density matrix and orbital analysis tools
* Integral screening and performance optimizations
* Extended basis set support
* Support for higher-order post-HF methods

#### Requirements
To build and run Hartree-Fock, you will need:

* A C++ compiler with C++23 support (e.g., GCC ≥ 13, Clang ≥ 16)
* CMake ≥ 3.26
* OpenMP (usually bundled with the compiler)
* A Unix-like environment (Linux or macOS or WSL)

<p align="justify" style="font-style: italic"> Note: The code has only been tested in Unix-like environments. Your mileage may vary if you are using Windows systems. If you are on a Windows Machine, it is highly recommended to install the package under WSL.</p>

#### Installation Instructions
Clone the repository and build using CMake:
```bash
git clone https://github.com/HemanthHaridas/hartee-fock-v2.git hartree-fock
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$(pwd)
cmake --build build 
cmake --install build 
```

To update the code:
```bash
cd hartree-fock
git pull
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$(pwd)
cmake --build build 
cmake --install build 
```

#### Usage Instructions

<p align="justify"> To run a Hartree–Fock calculation, provide an input file and redirect the output as desired: </p>

```bash 
hartee-fock input_file > output_file
``` 
<p align="justify"> All runtime information, including energies and convergence details, is written to standard output.

#### Input File
<p align="justify"> Planck uses a minimal, block-based input format inspired by traditional quantum chemistry codes. An example input file for single-point energy calculation on water at sto-3g basis set. </p>

```ini
[CALC]
BASIS       STO-3G
CALC_TYPE   ENERGY
THEORY      RHF
CHARGE      0 
MULTI       1
[END CALC]

[GEOM]
3
O   0.000000   0.000000   0.000000
H   0.757160   0.586260   0.000000
H  -0.757160   0.586260   0.000000
[END GEOM]
```

#### Geomerty Block

* The first line specifies the number of atoms
* Each subsequent line contains: ``` Element x y z```
* Coordinates are assumed to be: **Cartesian** and **Angstroms**

#### Calculation Block

| Keyword     | Description                          | Default Values |
|:-----------:|:------------------------------------:|:--------------:|
| `BASIS`     | Basis set name (e.g., `STO-3G`)      | `STO-3G`       |
| `CALC_TYPE` | Type of calculation (`ENERGY`, etc.) | `ENERGY`       |
| `THEORY`    | Electronic structure method (`RHF`)  | `RHF`          |
| `CHARGE`    | Total molecular charge               | `0`            |
| `MULTI`     | Spin multiplicity (2S + 1)           | `1`            |
| `USE_SYMM`  | Use point-group symmetry             | `ON`           |
| `USE_DIIS`  | Use DIIS in SCF cycles               | `ON`           |
| `MAXSCF`    | Maximum Number of SCF cycles         | `100`          |
| `TOLSCF`    | SCF Tolerance                        | `1E-10`        |
| `TOLERI`    | ERI Tolerance for Integral Screening | `1E-10`        |

#### Contributing to Hartree-Fock

<p align="justify"> Contributions to Hartree-Fock are welcome and encouraged, particularly those that improve clarity, correctness, performance, or educational value. This project is designed to be readable, hackable, and scientifically transparent, making it suitable for both learning and exploratory research in ab initio quantum chemistry.</p>


#### Ways to Contribute

You can contribute in many ways, including (but not limited to):

- Bug fixes and numerical correctness improvements  
- Performance optimizations (e.g., OpenMP parallelism, algorithmic improvements)  
- New features (e.g., UHF support, additional basis sets, analysis tools)  
- Documentation improvements (README, comments, theory notes)  
- Test cases and validation against reference data  
- Refactoring for improved structure and maintainability  

<p align="justify">Small contributions—such as improving comments or simplifying logic—are just as valuable as large feature additions.</p>

#### Development Guidelines

<p align="justify"> Please follow these guidelines when contributing to the codebase. </p>

#### Language and Standard

- Use **modern C++ (C++23)** features where appropriate  
- Avoid compiler-specific extensions unless absolutely necessary  

#### Code Style

- Prefer **readability over cleverness**  
- Use clear, descriptive variable and function names  
- Keep functions small and focused  
- Add comments where the physics, mathematics, or algorithms are non-trivial  

#### Scientific Correctness

- Clearly document assumptions and approximations  
- Where possible, reference standard literature (e.g., *Szabo & Ostlund*, *Helgaker et al.*)  
- Ensure numerical changes are validated against known results  

#### Parallelism

- OpenMP is used in performance-critical sections  
- Ensure thread safety and avoid hidden data races  
- Clearly document any parallel regions that rely on specific assumptions  

#### How to Contribute

1. Fork the repository on GitHub  

2. Create a new branch for your changes:
   ```bash
   git checkout -b feature/my-feature
   ```

3. Make your changes and commit them with clear, descriptive messages:
    ```bash
    git commit -m "Improve SCF convergence diagnostics"
    ```

4. Push your branch to your fork and open a Pull Request

When submitting a pull request, please include:

* A clear description of what the change does
* The motivation for the change
* Any relevant benchmarks, equations, or references

#### Reporting Issues
<p align="justify"> If you encounter a bug or unexpected behavior, please open an issue and include: </p>

- A minimal input example
- Relevant output or error messages
- Compiler version and platform details

<p align="justify"> Well-documented issues significantly reduce debugging time and are greatly appreciated. </p>

#### Project Philosophy
The project prioritizes:

- Correctness over features
- Clarity over opacity
- Educational value alongside performance

<p align="justify"> If you are unsure whether a proposed change aligns with the project’s direction, feel free to open an issue or a draft pull request for discussion.</p>



