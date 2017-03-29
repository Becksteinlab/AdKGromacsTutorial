.. -*- encoding: utf-8 -*-

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`


======================
Directory organization
======================

The workflow for setting up, running, and analysing a simulation
consists of multiple and rather different steps. It is useful to
perform these different steps in separate directories in order to
avoid overwriting files or using wrong files.

**Directory layout**

  For this tutorial the suggested directory layout is the following::

     coord/
     top/
     solvation/
     emin/
     posres/
     MD/
     analysis/

  You will work through these directories in sequence.

.. rubric:: Short description of the directories

:file:`coord`
  original PDB (structural) files
:file:`top`
  generating topology files (``.top``, ``.itp``)
:file:`solvation`
  adding solvent and ions to the system
:file:`emin`
  performing energy minimization
:file:`posres`
  short MD simulation with position restraints on the heavy protein
  atoms, to allow the solvent to equilibrate around the protein
  without disturbing the protein structure
:file:`MD`
  MD simulation (typically, you will transfer the :file:`md.tpr` file to a
  supercomputer, run the simulation there, then copy the the output
  back to this trajctory)
:file:`analysis`
  post-processing a production trajectory to facilitate easy visualization
  (i.e., using VMD); analysis of the simulations can be placed in
  (sub)directories under analysis, e.g. ::

     analysis/RMSD
     analysis/RMSF
     ...

  The subdirectories depend on the specific analysis tasks that you
  want to carry out. The above directory layout is only a suggestion
  but in practice some sort of ordered directory hierarchy has proven
  very useful.


.. Note ::

   The command snippets in this tutorial assume the above directory layout.
   The workflow is such that each step is carried out
   *inside the appropriate directory* and *relative* paths are used to
   access files from previous steps. It should be clear from the context
   in which directory the commands are to be executed. If you get a
   ``File input/output error`` from :program:`grompp` (or any of the
   other commands) then check that you are able to see the file by just
   doing a ``ls ../path/to/file`` from where you are in the file system.
   If you can't see the file then check (1) that you are in the correct
   directory, (2) that you have created the file in a previous step.


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
