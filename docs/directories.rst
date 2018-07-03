.. -*- encoding: utf-8 -*-

.. include:: /includes/defs.rst
.. include:: /includes/links.rst
   

======================
Directory organization
======================

The workflow for setting up, running, and analysing a simulation
consists of multiple and rather different steps. It is useful to
perform these different steps in separate directories in order to
avoid overwriting files or using wrong files.


Create working directories
==========================

It is recommended that the following directory structure be used, as the
tutorial steps through them sequentially::

     coord/
     top/
     solvation/
     emin/
     posres/
     MD/
     analysis/

Create these directories using::

     mkdir top solvation emin posres MD analysis


.. rubric:: Description of directories

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
  want to carry out. The above directory layout is only a suggestion,
  but, in practice, some sort of ordered directory hierarchy will facilitate
  reproducibility, improve efficiency, and maintain your sanity.

.. Important::

    The command snippets in this tutorial assume the directory layout given
    above as the workflow depends on each step's being carried out
    *inside the appropriate directory*. In particular, *relative* paths are used
    to access files from previous steps. It should be clear from context
    in which directory the commands are to be executed. If you get a
    ``File input/output error`` from :program:`grompp` (or any of the
    other commands), first check that you are able to see the file by just
    doing a ``ls ../path/to/file`` from where you are in the file system.
    If you can't see the file then check (1) that you are in the correct
    directory, (2) that you have created the file in a previous step.


Obtain starting structure
=========================

.. Note:: The starting structure :file:`coord/4ake_a.pdb` has been
          provided as part of the tutorial package, so the instructions that
          follow are optional for this tutorial. However, these steps provide an
          idea of what may be required in obtaining a suitable starting
          structure for MD simulation.

1. Download 4AKE_ the Protein Data Bank (PDB) through the web interface
2. Create a new PDB file with just chain A

   Modify the downloaded PDB file. For a relatively simple
   protein like AdK, one can just open the PDB file in a text editor and remove
   all the lines that are not needed.(For more complex situations, molecular
   modeling software can be used.)

  - Remove all comment lines (but keep TITLE, HEADER)
  - Remove all crystal waters (HOH) [#crystalwaters]_
  - Remove all chain B ATOM records.
  - Save as :file:`coord/4ake_a.pdb`.

.. rubric:: Footnotes

.. [#crystalwaters] Often you would actually want to retain
   crystallographic water molecules as they might have biological
   relevance. In our example this is likely not the case and by
   removing all of them we simplify the preparation step somewhat. If
   you keep them, :program:`gmx pdb2gmx` in the next step will
   actually create entries in the topology for them.


