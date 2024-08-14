# Invitation to the ESPResSo Summer School 2024

[![CECAM Flagship School registration link](https://img.shields.io/badge/CECAM%20Flagship%20School-Register%20Now-blue?style=for-the-badge)](https://www.cecam.org/workshop-details/1324)

The summer school "Simulating soft matter across scales" will take place on October 7-11, 2024, in Stuttgart. Registration is now open on [CECAM](https://www.cecam.org/workshop-details/1324).

# ESPResSo

[![GitLab CI](https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/badges/python/pipeline.svg)](https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/-/commits/python)
[![codecov](https://codecov.io/gh/espressomd/espresso/branch/python/graph/badge.svg)](https://codecov.io/gh/espressomd/espresso)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jngrad/espresso-binder/HEAD)
[![Contribute with Gitpod](https://img.shields.io/badge/Contribute%20with-Gitpod-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/espressomd/espresso)

This is the Molecular Dynamics software ESPResSo ("Extensible
Simulation Package for Research on Soft Matter Systems").

ESPResSo is an open-source platform designed for particle-based simulations using molecular dynamics, Monte Carlo techniques and couplings to lattice-based methods such as the lattice-Boltzmann method.
Focused on coarse-grained models, it can be used for research in soft matter science, statistical and biophysics, and process engineering.
To this end, a wide range of features is provided including, e.g., solvers for electrostatics and magnetostatics in various geometries, modelling chemical reactions and grand-canonical systems, hydrodynamics via Stokesian dynamics and the lattice-Boltzmann method, diffusion-advection-reaction equations, active particles and shear boundary conditions, and rigid body mechanics.

ESPResSo is controlled by a Python interface, allowing great flexibility in building and running simulation models. The simulation core is written in C++ and parallelized using MPI to allow for better performance.
Hence, the software can b employed on desktop machines, clusters as well as
on supercomputers. 

ESPResSo is used in scientific working groups all over the world both
as a production platform as well as a research platform for developing
new algorithms and methods and designing new models for coarse-grained
simulations.  It is mainly developed at the Institute for
Computational Physics of the University of Stuttgart, but has
contributors from all over the world.

## Documentation

* [User guide](https://espressomd.github.io/doc/index.html) Installation, use, and features
* [Tutorials](https://espressomd.github.io/#tutorials) for the simulation of various systems, e.g., Lennard-Jones particles, polymers and polyelectrolytes, particles in a fluid flow, ferrofluids, etc.
* [Video lectures](https://www.youtube.com/@espressosimulationpackage7832)
* [Documentation for developers](https://espressomd.github.io/#development)
* [Wiki](https://github.com/espressomd/espresso/wiki/)
* [Official website](https:/espressomd.org)

## Installation

You can try ESPResSo in the [cloud](https://espressomd.github.io/tutorials.html) without installing it. However, for good performance, install the software on your system.
 
Detailed installation instructions for Ubuntu and macOS can be found in the
user guide, section [Installation](https://espressomd.github.io/doc/installation.html).
Common installation issues are addressed in the
[FAQ](https://github.com/espressomd/espresso/wiki/Installation-FAQ).

For most users, we recommend downloading the latest release version of ESPResSo. You
can find it in the [release page](https://github.com/espressomd/espresso/releases),
together with past releases until 4.0. When choosing a release, we recommend that
you get the latest bugfix release in that line. For example, for 4.2 you would like
to use 4.2.2.

## Join the community

Please consider subscribing to our
[mailing list](https://espressomd.org/wordpress/community-and-support/mailing-lists/)
if you're actively using ESPResSo, as we occasionally need community
feedback when making decisions on the future of specific features in
upcoming releases. You'll also get notifications on bugfix releases.

## Please cite us!

If you use ESPResSo to publish scientific results, we would ask you to
acknowledge this usage by mentioning the software with its version number and
[citing the relevant papers](https://espressomd.org/wordpress/about/please-cite-us/).
A number of algorithms in ESPResSo are fairly advanced and unique to ESPResSo.
The authors of these contributions kindly ask you to cite the relevant
publications, as indicated in the documentation. For detailed instructions, see
[How to cite ESPResSo](https://espressomd.github.io/doc/introduction.html#how-to-cite-espresso).

## License

Copyright (C) 2010-2024 The ESPResSo project

Copyright (C) 2002-2010 Max-Planck-Institute for Polymer Research, Theory Group

ESPResSo is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

You should have received a [copy](COPYING) of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
