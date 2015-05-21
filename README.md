# QUMVIA
QUMVIA is a quantum mechanical software that implements Vibrational self-consistent field (VSCF) and vibrational configuration interaction (VCI) vibrational structure methods. VSCF and VCI are analogous to HF and CI methods in electronic structure theory, although for solving the vibrational wavefunction instead of the electronic one. 
The program need an external software for building a potential energy surface (PES) in the quartic force field representation (QFF). QUMVIA is able to use Gaussian, Gamess and Lio for this.

REQUIREMENTS
------------
* Intel MKL (Math Kernel Library).
* Intel Fortran compiler (can be obtained with a non-comercial licence).
* GNU Make.

COMPILATION
-----------
QUMVIA can be compilated in CPU version, which only permits the use of Gaussian or GAMESS for the PES buildup, or with in a Lio-compatible version, using GNU make. For example, the following compiles the CPU version:
```
make cpu
```
Available options:
* _cpu_: compiles QUMVIA version without Lio interface.
* _lio_: compiles the Lio-compatible version of the code.


