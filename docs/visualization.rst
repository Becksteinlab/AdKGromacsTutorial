.. -*- encoding: utf-8 -*-

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`


.. _trajectory-visualization:

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

We will use the `gmx trjconv`_ tool in Gromacs to center and remap our system.

.. Tip:: :program:`gmx trjconv` prompts the user with a number of questions that
         depend on the selected options. In the command line snippets below, the
         user input is directly fed to the standard input of :program:`trjconv`
         with the :kbd:`printf TEXT | gmx trjconv` "pipe" construct. In order to
         better understand the command, run it interactively without the pipe
         construct and manually provide the required information.

Center (:code:`-center`) on the *Protein* and remap all the molecules
(:code:`-pbc mol`) of the whole *System*::

  printf "Protein\nSystem\n" | gmx trjconv -s md.tpr -f md.xtc -center -ur compact -pbc mol -o md_center.xtc


Pinning down a tumbling protein
===============================

It is often desirable to *RMS-fit* the protein on a reference structure
(such as the first frame in the trajectory) to remove overall translation
and rotation. In Gromacs, the `gmx trjconv`_ tool can also do more "trajectory
conversion tasks". After (1) centering and remapping the system, we want
to (2) RMS-fit (due to technical limitations in :program:`gmx trjconv` you
cannot do both at the same time).

RMS-fit (:code:`-fit rot+trans`) to the protein *backbone* atoms in
the initial frame (supplied in the TPR file) and write out the
whole *System*::

  printf "Backbone\nSystem\n" | gmx trjconv -s md.tpr -f md_center.xtc -fit rot+trans -o md_fit.xtc


Check our modified trajectory
=============================

Visualize in VMD_::

  vmd ../posres/posres.pdb md_fit.xtc


.. Note:: If you don't have a :program:`vmd` command available on the command
          line then launch VMD_, load the ``posres/posres.pdb`` file
          (:menuselection:`File --> New Molecule...`), highlight your molecule 1
          ("em.pdb") and load the ``posres/md_fit.xtc`` trajectory into your
          *molecule 1*, :menuselection:`File --> Load Data Into Molecule`. You
          should see that the first frame (from the energy minimization) looks
          as if the water is in a distorted box shape whereas all further frames
          show a roughly spherical unit cell (the `rhombic dodecahedron`_).

	  
.. _Gromacs: http://www.gromacs.org  
.. _Gromacs manual: http://manual.gromacs.org/
.. _Gromacs documentation: http://www.gromacs.org/Documentation

.. _gmx trjconv: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-trjconv.html
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _rhombic dodecahedron: http://mathworld.wolfram.com/RhombicDodecahedron.html
