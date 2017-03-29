.. -*- encoding: utf-8 -*-

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`


===================
Energy minimization
===================

In order to remove "clashes" (i.e. close overlaps of the LJ cores) we
perform an *energy minimization*: Instead of a MD simulation we use an
algorithm to change the coordinates in such a way as to reduce the
total potential energy.


Telling Gromacs what it will do
===============================

First, we will copy a file from the templates folder (provided in this
tutorial) that tells Gromacs MD program *how* to do energy minimization::

  cp ../templates/em.mdp .

.. Note:: The MDP file :file:`em.mdp` is provided in
          :file:`templates/em.mdp`: copy it to the :file:`emin/`
          directory. You should have a look at it and modify it according to
          your needs. The individual parameters are explained under `mdp
          options`_.

This ``.mdp`` file contains the settings that dictate the nature of the
simulation. For energy minimization, we will use the simple *steepest
descent* minimizer (``integrator = steep`` in ``em.mdp``, which runs in
parallel). Use grompp_ (the GROMacs PreProcessor) to generate the run
input file (TPR) from the run parameter file (MDP), coordinate file
(the solvated system with ions; PDB), and the topology (TOP)::

  cd ../emin
  gmx grompp -f em.mdp -c ../solvation/ionized.pdb -p ../top/4ake.top -o em.tpr


Performing an energy minimization run
=====================================

The energy minimization is performed with :program:`mdrun` but by
using the appropriate ``integrator`` option in the `Run control
options in the MDP file`_ it has been instructed to do a energy
minimization::

  gmx mdrun -v -s em.tpr -deffnm em -c em.pdb

Ideally, the maximum force *Fmax* (gradient of the potential) should
be < 1e+03 |kJ/mol/nm**2| (but typically anything below 1e+05
|kJ/mol/nm**2| works). See the screen output or the :file:`em.log` file for
this information.

.. Note:: If you want to minimize further, you can use the :file"`em.pdb`
          structure as an input for a second run with either the *conjugate
          gradients* (``integrator = cg``) or the Newton-like
          *Broyden-Fletcher-Goldfarb-Shanno* (``integrator = l-bfgs``)
          minimizer. For details see `Run control options in the MDP file`_.


.. _`AdKTutorial.tar.bz2`:
    http://becksteinlab.physics.asu.edu/pages/courses/2013/SimBioNano/13/AdKTutorial.tar.bz2
.. _4AKE: http://www.rcsb.org/pdb/explore.do?structureId=4ake
.. _pdb2gmx: http://manual.gromacs.org/current/online/pdb2gmx.html
.. _editconf: http://manual.gromacs.org/current/online/editconf.html
.. _genbox: http://manual.gromacs.org/current/online/genbox.html
.. _genion: http://manual.gromacs.org/current/online/genion.html
.. _trjconv: http://manual.gromacs.org/current/online/trjconv.html
.. _trjcat: http://manual.gromacs.org/current/online/trjcat.html
.. _eneconv: http://manual.gromacs.org/current/online/eneconv.html
.. _grompp: http://manual.gromacs.org/current/online/grompp.html
.. _mdrun: http://manual.gromacs.org/current/online/mdrun.html
.. _`mdp options`: http://manual.gromacs.org/current/online/mdp_opt.html
.. _`Run control options in the MDP file`: http://manual.gromacs.org/current/online/mdp_opt.html#run
.. _`make_ndx`: http://manual.gromacs.org/current/online/make_ndx.html
.. _`g_tune_pme`: http://manual.gromacs.org/current/online/g_tune_pme.html
.. _gmxcheck: http://manual.gromacs.org/current/online/gmxcheck.html

.. _Gromacs manual: http://manual.gromacs.org/
.. _Gromacs documentation: http://www.gromacs.org/Documentation
.. _`Gromacs 4.5.6 PDF`: http://www.gromacs.org/@api/deki/files/190/=manual-4.5.6.pdf
.. _manual section: http://www.gromacs.org/Documentation/Manual

.. _`g_rms`: http://manual.gromacs.org/current/online/g_rms.html
.. _`g_rmsf`: http://manual.gromacs.org/current/online/g_rmsf.html
.. _`g_gyrate`: http://manual.gromacs.org/current/online/g_gyrate.html
.. _`g_dist`: http://manual.gromacs.org/current/online/g_dist.html
.. _`g_mindist`: http://manual.gromacs.org/current/online/g_mindist.html
.. _`do_dssp`: http://manual.gromacs.org/current/online/do_dssp.html

.. _DSSP: http://swift.cmbi.ru.nl/gv/dssp/
.. _`ATOM record of a PDB file`: http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
