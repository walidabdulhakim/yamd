# Molecular Dynamics in C++

Repository for coursework concerning the summer 2024 course [Molecular Dynamics in C++](https://pastewka.github.io/MolecularDynamics/) at the University of Freiburg, forked from [the MD Meson Template](https://github.com/imtek-simulation/meson-skeleton/).


- `/src` and `/tests` contain the library source code and tests in C++ using Eigen and gtest
- `/milestones` contain the source code for the binaries used to generate the data in the report

# To Run

```bash
cd <your repository>

# Configure and create build directory
meson setup builddir

# Compile
cd builddir
meson compile

# Run desired milestone
./milestones/01/01 
```
To get plots please run the jupyter files in each milestones after the milestones are run.
Also Note that the flag
`--buildtype=debug` should be changed to
`--buildtype=release`

milestone 09 is currently setup for small whisker and runs automatically for 2 temperature 0 and 600 and strainrates of 1e8 and 5e8 and can be changed to what every desired settings needed all one must do is go to the file and change it.

and one can run milestone 09 using mpi as such

```bash
mpirun -n 2 ./milestones/09/09

```

<<<<<<< Updated upstream
Huge credit to https://github.com/arulsingh, https://github.com/JulianKarrer and https://github.com/Dwaipayan-R-C as they inspired some of my codes 
=======
Huge credit to https://github.com/arulsingh, https://github.com/JulianKarrer and https://github.com/Dwaipayan-R-C as they inspired some of my codes 
>>>>>>> Stashed changes
