.. -*- encoding: utf-8 -*-

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`

==============================
Equilibrium molecular dynamics
==============================

Set up the production run
=========================

As usual, we must tell Gromacs what it will be doing using `gmx grompp`_
before we can perform our production simulation. Since we want to
start our run where we left off (after doing equilibration), we
prepare the TPR input file based on the last frame of the
:ref:`position-restraints` with :program:`gmx grompp`::

  cd ../MD
  cp ../templates/md.mdp .
  gmx grompp -f md.mdp -p ../top/4ake.top -c ../posres/posres.pdb -o md.tpr -maxwarn 3

The :file:`md.mdp` file uses different algorithms from the
:ref:`position-restraints` for the temperature and pressure coupling,
which are known to reproduce the exact *NPT* ensemble distribution.


Run the simulation
==================

Run the simulation as usual with `gmx mdrun`_::

  gmx mdrun -v -stepout 10 -s md.tpr -deffnm md -cpi

This will automatically utilize all available cores and GPUs if
available. The :code:`-cpi` flag indicates that you want Gromacs to
continue from a previous run. You can kill the job with
:kbd:`CONTROL-C`, look at the output, then continue with exactly the
same command line ::

  gmx mdrun -v -stepout 10 -s md.tpr -deffnm md -cpi

(Try it out!). The :code:`-cpi` flag can be used on the first run
without harm. For a continuation to occur, Gromacs needs to find the
checkpoint file :file:`md.cpt` and all output files (:file:`md.xtc`,
:file:`md.edr`, :file:`md.log`) in the current directory.


.. _gmx grompp: http://manual.gromacs.org/programs/gmx-grompp.html
.. _gmx mdrun: http://manual.gromacs.org/programs/gmx-mdrun.html
