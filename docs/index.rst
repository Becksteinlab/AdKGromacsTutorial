.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`

.. αβγδΔ


====================
AdK Gromacs Tutorial
====================

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


Tutorial files
==============

All of the necessary tutorial files can be found on GitHub and can be obtained
by cloning the repository::

  git clone https://github.com/Becksteinlab/AdKGromacsTutorial.git


Workflow overview
=================

For this tutorial we'll use Gromacs_ (version 5.1.3) to set up the system, run
the simulation, and perform analysis. An initial structure is provided, which
can be found in the :file:`tutorial/templates` directory, as well as the MDP
files that are necessary for input to Gromacs. The overall workflow consists of
the following steps:




.. _Gromacs: http://www.gromacs.org
.. _PDB: http://www.rcsb.org/pdb/home/home.do
