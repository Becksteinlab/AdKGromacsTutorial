.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`

.. αβγδΔ


====================
AdK Gromacs Tutorial
====================

|docs|

.. image:: docs/figs/adk_secondary.*
   :width: 30%
   :alt: Adenylate Kinase (AdK)
   :align: right

..   Adenylate Kinase (AdK). Secondary structure elements are colored
..   (magenta: α-helices, yellow: β-sheets).

The full tutorial is available online through
`Read the Docs<http://adkgromacstutorial.readthedocs.io>`__.

Objective
=========

Perform an all-atom molecular dynamics (MD) simulation—using the Gromacs_
MD package—of the apo enzyme adenylate kinase (AdK) in its open conformation in
a physiologically realistic environment, and carry out a basic analysis of its
structural properties in equilibrium.


Tutorial files
==============

All of the necessary tutorial files can obtained by cloning the repository::

  git clone https://github.com/Becksteinlab/AdKGromacsTutorial.git


Workflow overview
=================

For this tutorial we'll use Gromacs_ (version 5.1.3) to set up the system, run
the simulation, and perform analysis. An initial structure is provided, which
can be found in the :file:`tutorial/templates` directory, as well as the MDP
files that are necessary for input to Gromacs. The overall workflow consists of
the following steps:

1. Download tutorial files and set up working directories

   - Obtain structure 4AKE from PDB_
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


.. _Gromacs: http://www.gromacs.org
.. _PDB: http://www.rcsb.org/pdb/home/home.do

.. |docs| image:: https://readthedocs.org/projects/adkgromacstutorial/badge/?version=master
    :alt: Documentation Status
    :scale: 100%
    :target: https://readthedocs.org/projects/adkgromacstutorial
