.. -*- encoding: utf-8 -*-

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`

========
Analysis
========

A typical analysis tasks reads the trajectory (XTC) or energy (EDR)
file, computes quantities, and produces datafiles that can be plotted
or processed further, e.g. using Python scripts. A strength of
Gromacs_ is that it comes with a wide range of tools that each do one
particular analysis task well (see the `Gromacs manual`_ and the
`Gromacs documentation`_).


Observables
-----------

A number of interesting quantities and observables [#observable]_ can
be calculated with Gromacs tools. A selection is shown below but you
are encouraged to read the `Gromacs manual`_ and the `Gromacs
documentation`_ to find out what else is available.

.. rubric:: Selection of Gromacs analysis tools

`gmx energy <g_energy>`_
    basic thermodynamic properties of the system

`gmx rms <g_rms>`_
    calculate the root mean square deviation from a reference structure

`gmx rmsf <g_rmsf>`_
    calculate the per-residue root mean square fluctuations

`gmx gyrate <g_gyrate>`_
    calculate the radius of gyration

`gmx mindist <g_mindist>`_
    calculate the distance between atoms or groups of atoms (make a
    index file with `make_ndx`_ to define the groups of
    interest). :program:`gmx mindist` is especially useful to find water
    molecules close to a region of interest.

    For AdK, look at the distance between the |Calpha| of K145 and I52.

`gmx do_dssp <do_dssp>`_
    Use the DSSP_ algorithm [Kabsch1983]_ to analyze the secondary structure
    (helices, sheets, ...). 


Basic analysis
--------------

.. toctree::
   :maxdepth: 1
   :numbered:

   analysis/rmsd
   analysis/rmsf
   analysis/distances
   analysis/rgyr
   
    
.. See also::

   * For analysis coupled with visualization look at VMD_.
   * For analyzing MD trajectories in many common formats (including
     the XTC, TRR, etc. used by Gromacs) using Python_, have a look at
     the `MDAnalysis`_ Python library (the `MDAnalysis Tutorial`_ is a
     good place to start... and it also uses AdK as an example).


.. rubric:: Footnotes

.. [#observable] "Observable" is used in the widest sense in that we
   know an estimator function of all or a subset of the system's phase
   space coordinates that is averaged to provide a quantity of
   interest. In many cases it requires considerable more work to
   connect such an "observable" to a true experimental observable that
   is measured in an experiment.

     

.. _MDAnalysis: http://www.mdanalysis.org/
.. _`MDAnalysis Tutorial`: http://www.mdanalysis.org/MDAnalysisTutorial/
.. _Python: https://www.python.org

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

.. _Gromacs: http://www.gromacs.org
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
