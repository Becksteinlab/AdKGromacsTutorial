.. -*- encoding: utf-8 -*-

.. include:: /includes/defs.rst
.. include:: /includes/links.rst


==================================
Generate a solvated protein system
==================================


Generate a topology
===================

Using the modified PDB file (chain A of 4AKE with crystal waters removed),
generate a topology file for the CHARMM27 force field together with the
TIP3P water model using the :ref:`gmx pdb2gmx` tool::

    cd top
    gmx pdb2gmx -f ../coord/4ake_a.pdb -o protein.pdb -p 4ake.top -i protein_posre.itp -water tip3p -ff charmm27

.. Note:: Total charge -4.000e (in the next step we will add ions to
          neutralize the system; we need a net-neutral system to properly
          handle electrostatics)


Solvate the protein
===================

Adding water
------------

Create a simulation box with :ref:`gmx editconf` and add solvent with
:ref:`gmx solvate`::

  cd ../solvation
  gmx editconf -f ../top/protein.pdb -o boxed.pdb -c -d 0.5 -bt dodecahedron
  gmx solvate -cp boxed.pdb -cs spc216 -p ../top/4ake.top -o solvated.pdb

.. Attention:: In order to reduce the system size and make the simulations run
               faster we are choosing a very tight box (minimum protein-edge
               distance 0.5 nm, ``-d 0.5``); for simulations you want to publish
               this number should be 1.2...1.5 nm so that the electrostatic
               interactions between copies of the protein across periodic
               boundaries are sufficiently screened.

:program:`gmx solvate` updates the number of solvent molecules ("SOL") in the
topology file (check the ``[ system ]`` section in
``top/system.top``) [#topupdate]_.


Adding ions
-----------

Ions can be added with the :ref:`gmx genion` program in Gromacs.

First, we need a basic TPR file (an empty file is sufficient, just
ignore the warnings that :ref:`gmx grompp` spits out by setting
``-maxwarn 10``), then run :ref:`gmx genion` (which has convenient
options to neutralize the system and set the concentration (check
the help!); :ref:`gmx genion` also updates the topology's ``[ system
]`` section if the top file is provided [#topupdate]_; it reduces the
"SOL" molecules by the number of removed molecules and adds the
ions, e.g. "NA" and "CL"). ::

  touch ions.mdp
  gmx grompp -f ions.mdp -p ../top/4ake.top -c solvated.pdb -o ions.tpr
  printf "SOL" | gmx genion -s ions.tpr -p ../top/4ake.top -pname NA -nname CL -neutral -conc 0.15 -o ionized.pdb

The final output is :file:`solvation/ionized.pdb`. Check visually in
VMD_ (but note that the dodecahedral box is not represented
properly). [#visualization]_.


.. rubric:: Footnotes

.. [#topupdate] The automatic modification of the top file by
   :ref:`gmx solvate` and :ref:`gmx genion` can become a
   problem if you try to run these commands multiple times and you get
   error messages later (typically from :ref:`gmx grompp`) that
   the number of molecules in structure file and the topology file do
   not agree. In this case you might have to manually delete or adjust
   the corresponding lines in :file:`system.top` file.

.. [#visualization] For notes on how to visualize MD systems with VMD_
   see the notes on :ref:`trajectory-visualization`.
