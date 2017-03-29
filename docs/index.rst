.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`

.. αβγδΔ


==================
AdKGromacsTutorial
==================

.. image:: /figs/adk_secondary.*
   :width: 30%
   :alt: Adenylate Kinase (AdK)
   :align: right

..   Adenylate Kinase (AdK). Secondary structure elements are colored
..   (magenta: α-helices, yellow: β-sheets).


Objective
=========

Perform an all-atom molecular dynamics (MD) simulation—using the Gromacs_
MD package—of the apo enzyme adenylate kinase (AdK) in its open conformation in
a physiologically realistic environment, and carry out a basic analysis of its
structural properties in equilibrium.


Workflow overview
=================

For this tutorial we'll use Gromacs_ (version 5.1.3) to set up the system, run
the simulation, and perform analysis. The overall workflow consists of the
following steps:

1. Download tutorial files and set up working directories

2. Prepare the system

   - Obtain structure 4AKE from PDB_, select chain A

   - Generate topology using default protonation states

   - Solvate in water in simulation cell (rhombic dodecahedron)

   - Add NaCl ions to neutralize and obtain final physiological concentration
     of 100 mM

3. Perform energy minimization (:file:`tutorial/emin`)

4. Perform position-restrained equilibration in the in the isobaric-isothermal
   (NPT) ensemble (:file:`tutorial/posres`)

5. Run a production MD simulation in the NPT ensemble (:file:`tutorial/MD`)

6. Visualize the trajectory

   - Center the protein in the box (periodic boundary conditions)

   - RMS-fit the protein in each snapshot to the first snapshot

This tutorial provides an initial structure, which can be found in the
:file:`tutorial/coord` directory, as well as the necessary MDP files that are
necessary for input to Gromacs in :file:`tutorial/templates`.


.. _Gromacs: http://www.gromacs.org
.. _PDB: http://www.rcsb.org/pdb/home/home.do


--------------------------------------------------------------------------------

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Contents

   directories
   preparation
   emin
   posres
   simulation
   visualization
   analysis


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
