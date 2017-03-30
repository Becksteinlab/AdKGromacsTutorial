.. -*- encoding: utf-8 -*-

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`

==============================
Equilibrium molecular dynamics
==============================

Set up the production run
=========================

As usual, we must tell Gromacs what it will be doing using grompp_
before we can perform our production simulation. Since we want to
start our run where we left off (after doing equilibration), we
prepare the TPR input file based on the last frame of the
:ref:`position-restraints` with grompp_::

  cd ../MD
  cp ../templates/md.mdp .
  gmx grompp -f md.mdp -p ../top/4ake.top -c ../posres/posres.pdb -o md.tpr -maxwarn 3

The :file:`md.mdp` file uses different algorithms from the
:ref:`position-restraints` for the temperature and pressure coupling,
which are known to reproduce the exact *NPT* ensemble distribution.


Run the simulation
==================

.. rubric:: CPU run

If your workstation has a decent number of cores or if you simply
don't mind waiting a bit longer you can also run the simulation as
usual with ::

  gmx mdrun -v -stepout 10 -s md.tpr -deffnm md -cpi

This will automatically utilize all available cores. The :code:`-cpi`
flag indicates that you want Gromacs to continue from a previous
run. You can kill the job with :kbd:`CONTROL-C`, look at the output,
then continue with exactly the same command line ::

  gmx mdrun -v -stepout 10 -s md.tpr -deffnm md -cpi

(Try it out!). The :code:`-cpi` flag can be used on the first run
without harm. For a continuation to occur, Gromacs needs to find the
checkpoint file :file:`md.cpt` and all output files (:file:`md.xtc`,
:file:`md.edr`, :file:`md.log`) in the current directory.


.. rubric:: GPU run

We can also try utilizing the GPU(s) available on the workstation. We use
a modified MDP file that contains settings compatible with GPU-based
Gromacs simulations to generate a new TPR file, which is used to perform
the simulation::

  gmx grompp -f mdgpu.mdp -p ../top/4ake.top -c ../posres/posres.pdb -o mdgpu.tpr -maxwarn 3
  gmx mdrun -v -stepout 10 -s mdgpu.tpr -deffnm mdgpu -cpi

On machines equipped with a well-matched CPU and GPU, a GPU-accelerated Gromacs
run can be around 3–5 times faster than a CPU-only run on the same machine.


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

.. _saguaro: http://a2c2.asu.edu/resources/saguaro/
.. _`How to login to saguaro`: http://a2c2.asu.edu/how-to-2/
.. _ASU: http://asu.edu
.. _`PHY494/PHY598/CHM598 — Simulation approaches to Bio- and Nanophysics`:
   http://becksteinlab.physics.asu.edu/learning/28/phy494-phy598-chm598-simulation-approaches-to-bio-and-nanophysics
