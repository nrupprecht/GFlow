# GFlow v 4.0

## Getting Started 

Fourth generation of the granular flow simulator.

This generation was heavily inspired by LAMMPS, GROMACS, OpenMD, etc. to make improvements to the program structure. For example, most objects inherit from a base class that has pointers to the GFlow simulation object and all of its main objects (inspired by LAMMPS).

This version also includes Doxygen style commenting which can be used to generate a users manual or reference.

## Useful sites

[LAMMPS](https://github.com/lammps/lammps)

## Resources

* Pall S, Hess B. Comp Phys Comm 2013;184:2641–50. - Interesting new take on Verlet lists. Used in GROMACS.

* Vectorized forces/neighbor lists [here](ftp://crack.seismo.unr.edu/downloads/russell/O(N)/grest_1989_vectorized_link_cell_code_md.PDF)

* Parallel decompositions [paper](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.35.6047&rep=rep1&type=pdf)

* GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers - Paper on GROMACS

* ls1 mardyne: Current record holder for largest md simulation - over 10^12 particles. See this [paper](https://arxiv.org/pdf/1408.4599.pdf)

* OpenMD: See e.g. [this](http://openmd.org/wp-content/docs/OpenMD-2.5.pdf)

Run with
```
./bin/driver
```

## TODO:

* Bodies, rigid and otherwise.

* Walls

* Angles / torque

* Sectorization that is optimized for multiple scales of particles

* Hard triangles

* Bond/angle modifier

* Time step monitoring --- Completed.

* Vectorize integration step

* Adding/removing particles while the simulation is running

* More efficient force calculation

* Parallelization

## Authors
* **Nathaniel Rupprecht** - *Author*