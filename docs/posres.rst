.. -*- encoding: utf-8 -*-

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`

.. _position-restraints:

=================================
Position-restrained equilibration
=================================

We first perform a short MD simulation with harmonic position
restraints on the heavy protein atoms. This allows the solvent to
equilibrate around the protein without disturbing the protein
structure. In addition, we use "weak coupling" temperature and
pressure coupling algorithms to obtain the desired temperatue,
:math:`T = 300` K, and pressure, :math:`P = 1` bar.


Set up and generate the run file
================================

We must first tell Gromacs *how* to perform our equilibration run
in the same way that we did for the energy minimization step.
This step requires the :file:`top/protein_posres.itp` file with the
default value for the harmonic force constants of 1000
|kJ/mol/nm**2|. The position restraints are switched on by setting the
:code:`-DPOSRES` flag in the :file:`posres.mdp` file (see `mdp options`_).

Create the run input (TPR) file, using the :ref:`energy minimized
system <energy-minimization>` as the starting structure::

  cd ../posres
  cp ../templates/posres.mdp .
  gmx grompp -f posres.mdp -o posres.tpr -p ../top/4ake.top -c ../emin/em.pdb -maxwarn 2

The mdp file contains cut-off settings that approximate the native
CHARMM values (in the CHARMM program).

Weak (Berendsen) coupling is used for both temperature and pressure to
quickly equilibrate. The protein and the solvent (water and ions) are
coupled as separate groups. Gromacs provides a range of groups
automatically (run :samp:`gmx make_ndx -f {TPR}` to see them) and we use
the groups ``Protein`` and ``non-Protein`` (these particularly groups
work since roughly Gromacs 4.5.3). If the standard groups do not work
then you will have to create the groups yourself using :samp:`gmx make_ndx
-f {TPR} -o md.ndx` (which would save them in a file :file:`md.ndx`) and
supply it to :code:`gmx grompp -n md.ndx`.


Perform equilibration
=====================

Run the position restraints equilibration simulation::

  gmx mdrun -v -stepout 10 -s posres.tpr -deffnm posres -c posres.pdb

.. Attention:: Here the runtime of 10 ps is too short for real production
               use; typically 1 to 5 ns are used.

	       
.. rubric:: Generate a centered trajectory in the primary unitcell

In order to visually check your system, first create trajectory with all
molecules in the primary unitcell (``-ur compact``; see also below the
more extensive notes on :ref:`trajectory-visualization`)::

   echo "System" | gmx trjconv -ur compact -s posres.tpr -f posres.xtc -pbc mol -o posres_ur.xtc

.. rubric:: Visually check centered trajectory in VMD

If you have VMD_ installed then you can quickly visualize the system
with the command ::	    
	    
   vmd ../emin/em.pdb posres_ur.xtc

   
If you don't have a :program:`vmd` command available on the command
line then launch VMD_, load the ``emin/em.pdb`` file
(:menuselection:`File --> New Molecule...`), highlight your molecule 1
("em.pdb") and load the ``posres/posres_ur.xtc`` trajectory into your
*molecule 1*, :menuselection:`File --> Load Data Into Molecule`. You
should see that the first frame (from the energy minimization) looks
as if the water is in a distorted box shape whereas all further frames
show a roughly spherical unit cell (the `rhombic dodecahedron`_).


.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _rhombic dodecahedron: http://mathworld.wolfram.com/RhombicDodecahedron.html

.. _`mdp options`: http://manual.gromacs.org/online/mdp_opt.html
