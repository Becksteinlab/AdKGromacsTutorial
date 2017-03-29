.. -*- encoding: utf-8 -*-

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`

========================
Trajectory visualization
========================

Analysis are normally performed locally on a workstation,
i.e. copy back all the files from the supercomputer to your local
directory.

A typical analysis tasks reads the trajectory (XTC) or energy (EDR)
file, computes quantities, and produces data files that can be plotted
or processed further, e.g. using Python scripts. A strength of
Gromacs_ is that it comes with a wide range of tools that each do one
particular analysis task well (see the `Gromacs manual`_ and the
`Gromacs documentation`_).


Keeping the protein in one piece
================================

If you just look at the output trajectory :file:`md.xtc` in VMD_ then
you will see that the protein can be split across the periodic
boundaries and that the simulation cell just looks like a distorted
prism. You should *recenter* the trajectory so that the protein is at
the center, *remap* the water molecules (and ions) to be located in a
more convenient unitcell representation.

We will use the trjconv_ tool in Gromacs to center and remap our system.

.. Note::

   :program:`trjconv` asks the user a number of
   questions that depend on the chosen options. In the command line
   snippets below, this user input is directly fed to the standard input
   of :program:`trjconv` with the :kbd:`printf TEXT | trjconv` "pipe"
   construct. In order to better understand the command, run it
   interactively without the pipe construct and manually provide the
   required information.

Center (:code:`-center`) on the *Protein* and remap all the molecules
(:code:`-pbc mol`) of the whole *System*::

  printf "Protein\nSystem\n" | gmx trjconv -s md.tpr -f md.xtc -center -ur compact -pbc mol -o md_center.xtc


Pinning down a tumbling protein
===============================

It is often desirable to *RMS-fit* the protein on a reference structure
(such as the first frame in the trajectory) to remove overall translation
and rotation. In Gromacs, the trjconv_ tool can also do more "trajectory
conversion tasks". After (1) centering and remapping the system, we want
to (2) RMS-fit (due to technical limitations in :program:`trjconv` you
cannot do both at the same time).

RMS-fit (:code:`-fit rot+trans`) to the protein *backbone* atoms in
the initial frame (supplied in the TPR file) and write out the
whole *System*::

  printf "Backbone\nSystem\n" | gmx trjconv -s md.tpr -f md_center.xtc -fit rot+trans -o md_fit.xtc


Check our modified trajectory
=============================

Visualize in VMD_::

  vmd ../posres/posres.pdb md_fit.xtc



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

.. _`g_energy`: http://manual.gromacs.org/current/online/g_energy.html
.. _`g_rms`: http://manual.gromacs.org/current/online/g_rms.html
.. _`g_rmsf`: http://manual.gromacs.org/current/online/g_rmsf.html
.. _`g_gyrate`: http://manual.gromacs.org/current/online/g_gyrate.html
.. _`g_dist`: http://manual.gromacs.org/current/online/g_dist.html
.. _`g_mindist`: http://manual.gromacs.org/current/online/g_mindist.html
.. _`do_dssp`: http://manual.gromacs.org/current/online/do_dssp.html

.. _DSSP: http://swift.cmbi.ru.nl/gv/dssp/
.. _`ATOM record of a PDB file`: http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

.. _saguaro: http://a2c2.asu.edu/resources/saguaro/
.. _Gromacs: http://www.gromacs.org
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
