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

`gmx energy`_
    basic thermodynamic properties of the system

`gmx rms`_
    calculate the root mean square deviation from a reference structure

`gmx rmsf`_
    calculate the per-residue root mean square fluctuations

`gmx gyrate`_
    calculate the radius of gyration

`gmx mindist`_
    calculate the distance between atoms or groups of atoms (make a
    index file with `gmx make_ndx`_ to define the groups of
    interest). :program:`gmx mindist` is especially useful to find water
    molecules close to a region of interest.

    For AdK, look at the distance between the |Calpha| of K145 and I52.

`gmx do_dssp`_
    Use the DSSP_ algorithm [Kabsch1983]_ to analyze the secondary structure
    (helices, sheets, ...). 


Basic analysis
--------------

.. toctree::
   :maxdepth: 1

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

     

.. updated links (2018)

.. Gromacs
.. _Gromacs: http://www.gromacs.org   

.. _Gromacs manual: http://manual.gromacs.org/
.. _Gromacs documentation: http://www.gromacs.org/Documentation
.. _manual section: http://www.gromacs.org/Documentation/Manual

.. _gmx pdb2gmx: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html
.. _gmx editconf: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-editconf.html
.. _gmx solvate: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-solvate.html
.. _gmx genion: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-genion.html
.. _gmx trjconv: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-trjconv.html
.. _gmx trjcat: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-trjcat.html
.. _gmx eneconv: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-eneconv.html
.. _gmx grompp: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-grompp.html
.. _gmx mdrun: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-mdrun.html
.. _`mdp options`: http://manual.gromacs.org/online/mdp_opt.html
.. _`Run control options in the MDP file`: http://manual.gromacs.org/online/mdp_opt.html#run
.. _`gmx make_ndx`: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-make_ndx.html
.. _`gmx tune_pme`: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-tune_pme.html
.. _gmx check: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-check.html
.. _`gmx energy`: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-energy.html
.. _`gmx rms`: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-rms.html
.. _`gmx rmsf`: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-rmsf.html
.. _`gmx gyrate`: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-gyrate.html
.. _`gmx distance`: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-dist.html
.. _`gmx mindist`: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-mindist.html
.. _`gmx do_dssp`: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-do_dssp.html


.. Protein Databank   
.. _4AKE: http://www.rcsb.org/pdb/explore.do?structureId=4ake

.. _`ATOM record of a PDB file`: http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

.. Other programs
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _DSSP: http://swift.cmbi.ru.nl/gv/dssp/
.. _Python: https://www.python.org

.. MDAnalysis (shameless plug)
.. _MDAnalysis: http://www.mdanalysis.org/
.. _`MDAnalysis Tutorial`: http://www.mdanalysis.org/MDAnalysisTutorial/

