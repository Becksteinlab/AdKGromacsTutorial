.. AdKGromacsTutorial documentation master file, created by
   sphinx-quickstart on Tue Mar 28 16:31:11 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`

.. αβγδΔ


AdKGromacsTutorial: Simulating AdK with Gromacs
===============================================

.. image:: /figs/adk_secondary.*
   :width: 30%
   :alt: Adenylate Kinase (AdK)
   :align: right

..   Adenylate Kinase (AdK). Secondary structure elements are colored
..   (magenta: α-helices, yellow: β-sheets).


Objective
---------
Perform aa MD simulation of the the enzyme adenylate kinase
(AdK) in its open conformation without a ligand bound. Simulate it in
a realistic environment (100 mM NaCl solution at :math:`T = 300` K and
:math:`P = 1` bar) and analyze its structural properties.

This tutorial progresses through the individual steps needed to set up
and run an equilibrium MD simulation of AdK using Gromacs.


Contents:

.. rubric:: Workflow summary
.. toctree::
   :maxdepth: 1
   :numbered:

   directory_organization
   system_setup
   energy_minimization
   position_restraints_MD
   equilibrium_MD
   trajectory_visualization


Overview of workflow
--------------------
For this tutorial we'll use Gromacs_ (version 5.1.3) to set up the system, run the simulation, and perform analysis. The overall work flow contains the following steps:

  1. Download tutorial files and organize the work space

  2. Setup

     - Obtain structure 4AKE from PDB, select chain A

     - Use default protonation states

     - Generate topology

     - Solvate in water in simulation cell (rhombic dodecahedron)

     - Add NaCl ions to neutralize and final physiological concentration
       of 100 mM

  3. Energy minimization (EM)

  4. Position restraint equilibration of solvent (MD); *NPT* (use weak
    coupling (Berendsen) schemes)

  5. Equilibrium MD simulation (unrestrained, *NPT*, use Nose-Hoover and
     Parrinello-Rahman temperature and pressure coupling)

  6. Trajectory visualization

     - Center the protein in the box (periodic boundary conditions)
     - RMS-fit the protein in each snapshot to the first snapshot

  All input files are provided in the same directory as :file:`AdKTutorial.html`. Start by uncompressing the package file::

    tar -jxvf AdKTutorial.tar.bz2
    cd AdKTutorial

  A starting structure can be found in the :file:`tutorial/coord` directory and MDP files are in :file:`AdKTutorial/templates`.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _Gromacs: http://www.gromacs.org
