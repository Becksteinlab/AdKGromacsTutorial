.. -*- encoding: utf-8 -*-

.. include:: /includes/defs.rst
.. include:: /includes/links.rst

.. _energy-minimization:		      

===================
Energy minimization
===================

In order to remove "clashes" (i.e. close overlaps of the LJ cores) we
perform an *energy minimization*: Instead of a MD simulation we use an
algorithm to change the coordinates in such a way as to reduce the
total potential energy.


Set up and generate the run file
================================

First, we will copy a file from the templates folder (provided in this
tutorial) that tells Gromacs MD program *how* to do energy minimization::

  cp ../templates/em.mdp .

.. Tip:: Have a look at the MDP file to get a feel for what kinds of settings
         can be adjusted to suit one's needs. Individual parameters are
         explained in more detail in `mdp options`_.

The ``*.mdp`` file contains the settings that dictate the nature of the
simulation. For energy minimization, we will use the simple *steepest
descent* minimizer (``integrator = steep`` in ``em.mdp``, which runs in
parallel). Use :ref:`gmx grompp` (the GROMacs PreProcessor) to generate the run
input file (TPR) from the run parameter file (MDP), coordinate file
(the solvated system with ions; PDB), and the topology (TOP)::

  cd ../emin
  gmx grompp -f em.mdp -c ../solvation/ionized.pdb -p ../top/4ake.top -o em.tpr


Perform energy minimization
===========================

The energy minimization is performed with :ref:`gmx mdrun` but by
using the appropriate ``integrator`` option in the `Run control
options in the MDP file`_ it has been instructed to do a energy
minimization::

  gmx mdrun -v -s em.tpr -deffnm em -c em.pdb

Ideally, the maximum force *Fmax* (gradient of the potential) should
be < 1e+03 |kJ/mol/nm**2| (but typically anything below 1e+05
|kJ/mol/nm**2| works). See the screen output or the :file:`em.log` file for
this information.

.. Tip:: The final frame of minimization (the structure in :file:`em.pdb`) can
         be used as the input structure for further minimization runs. It is
         common to do an initial energy minimization using the efficient
         steepest descent method and further minimization with a more
         sophisticated method such as the *conjugate gradient* algorithm
         (``integrator = cg``) or the Newton-like
         *Broyden-Fletcher-Goldfarb-Shanno* (``integrator = l-bfgs``) minimizer.
         For details, see `Run control options in the MDP file`_.

