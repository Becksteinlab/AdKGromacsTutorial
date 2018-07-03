.. -*- encoding: utf-8 -*-

.. include:: /includes/defs.rst
.. include:: /includes/links.rst
.. highlight:: bash
	     
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


:Gromacs: aiming for 2018 (maybe will work for 5.x, 2016)
:Tutorial: |release|
:Date: |today|

.. warning::

   This tutorial was originally based on an `older tutorial for
   Gromacs 4.x
   <https://becksteinlab.physics.asu.edu/learning/33/tutorial-simulating-adk-with-gromacs>`_
   and has not been completely transitioned to modern Gromacs
   versions. It will not work seamlessly and a number of MDP options
   are outdated and need to be updated. Please raise any problems in
   the `issue tracker
   <https://github.com/Becksteinlab/AdKGromacsTutorial/issues>`_.

   
.. seealso::

   Justin Lemkul's excellent `GROMACS Tutorials
   <http://www.mdtutorials.com/gmx/index.html>`_, which have recently
   been updated for Gromacs 2018.
       

Objective
=========

Perform an all-atom molecular dynamics (MD) simulation—using the Gromacs_
MD package—of the apo enzyme adenylate kinase (AdK) in its open conformation in
a physiologically realistic environment, and carry out a basic analysis of its
structural properties in equilibrium.


Tutorial files
==============

All of the necessary tutorial files can be found on GitHub in the
`Becksteinlab/AdKGromacsTutorial/tutorial
<https://github.com/Becksteinlab/AdKGromacsTutorial/tree/master/tutorial>`_
directory, which can be easily obtained by git-cloning the repository::

  git clone https://github.com/Becksteinlab/AdKGromacsTutorial.git


Workflow overview
=================

For this tutorial we'll use Gromacs_ (versions 5, 2016, 2018 should
work) to set up the system, run the simulation, and perform
analysis. An initial structure is provided, which can be found in the
:file:`tutorial/templates` directory, as well as the MDP files that
are necessary for input to Gromacs. The overall workflow consists of
the following steps:

.. toctree::
   :maxdepth: 1
   :numbered:

   directories
   preparation
   emin
   posres
   simulation
   visualization
   analysis
   references

.. _Gromacs: http://www.gromacs.org

