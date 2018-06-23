# GFlow v 4.0

## Getting Started 

Fourth generation of the granular flow simulator.

This generation was heavily inspired by LAMMPS, GROMACS, OpenMD, etc. to make improvements to the program structure. For example, most objects inherit from a base class that has pointers to the GFlow simulation object and all of its main objects (inspired by LAMMPS).

## Useful sites

[LAMMPS](https://github.com/lammps/lammps)

## Resources

* Pall S, Hess B. Comp Phys Comm 2013;184:2641–50. - Interesting new take on Verlet lists. Used in GROMACS.

* GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers - Paper on GROMACS

[29]

Run with
```
./bin/driver
```

## Authors
* **Nathaniel Rupprecht** - *Author*