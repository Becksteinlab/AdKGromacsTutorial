====================
AdK Gromacs Tutorial
====================

|docs|

![AdK](/docs/figs/adk_secondary.jpg?raw=true "Adenylate Kinase")

The full tutorial is available online through
[Read the Docs](http://adkgromacstutorial.readthedocs.io).

Objective
=========

Perform an all-atom molecular dynamics (MD) simulation—using the [Gromacs]
MD package—of the apo enzyme adenylate kinase (AdK) in its open conformation in
a physiologically realistic environment, and carry out a basic analysis of its
structural properties in equilibrium.


Tutorial files
==============

All of the necessary tutorial files can obtained by cloning the repository

```bash
  git clone https://github.com/Becksteinlab/AdKGromacsTutorial.git
```


Workflow overview
=================

For this tutorial we'll use [Gromacs] (version 5.1.3) to set up the system, run
the simulation, and perform analysis. An initial structure is provided, which
can be found in the :file:`tutorial/templates` directory, as well as the MDP
files that are necessary for input to Gromacs. The overall workflow consists of
the following steps:

1. Download tutorial files and set up working directories

   - Obtain structure 4AKE from [PDB]
   - Generate a stripped PDB file containing only chain A and no crystal waters

2. Solvate the protein system

   - Generate topology using default protonation states

   - Solvate in water in simulation cell (rhombic dodecahedron)

   - Ionize system with NaCl to neutralize and obtain physiological concentration

3. Perform energy minimization

4. Perform position-restrained equilibration

5. Run production MD in the NPT ensemble

6. Visualize the trajectory

   - Center the protein in a box with primitive unitcell representation (periodic boundary conditions)

   - RMS-fit the protein in each snapshot to the first snapshot


[Gromacs]: http://www.gromacs.org
[PDB]: http://www.rcsb.org/pdb/home/home.do

.. |docs| image:: https://readthedocs.org/projects/adkgromacstutorial/badge/
    :alt: Documentation Status
    :scale: 100%
    :target: https://readthedocs.org/projects/adkgromacstutorial
