.. -*- encoding: utf-8 -*-

.. include:: /includes/defs.rst
.. include:: /includes/links.rst
   

========
Analysis
========

A typical analysis tasks reads the trajectory (XTC) or energy (EDR)
file, computes quantities, and produces datafiles that can be plotted
or processed further, e.g. using Python scripts. A strength of
Gromacs_ is that it comes with a wide range of tools that each do one
particular analysis task well (see the `Gromacs manual`_ and the
`Gromacs documentation`_).


Basic analysis
--------------

As examples, we perform a number of common analysis tasks.

.. toctree::
   :maxdepth: 1

   analysis/rmsd
   analysis/rmsf
   analysis/distances
   analysis/rgyr


More Gromacs tools
------------------

A number of interesting quantities and observables [#observable]_ can
be calculated with Gromacs tools. A selection is shown below but you
are encouraged to read the `Gromacs manual`_ and the `Gromacs
documentation`_ to find out what else is available.
	      
.. rubric:: Selection of Gromacs analysis tools

The full `list of Gromacs commands`_ contains 98 different tools. A
small selection of commonly used ones are shown here:

:ref:`gmx energy`
    basic thermodynamic properties of the system

:ref:`gmx rms`
    calculate the root mean square deviation from a reference structure

:ref:`gmx rmsf`
    calculate the per-residue root mean square fluctuations

:ref:`gmx gyrate`
    calculate the radius of gyration

:ref:`gmx mindist` and :ref:`gmx distance`
    calculate the distance between atoms or groups of atoms (make a
    index file with :ref:`gmx make_ndx` to define the groups of
    interest). :program:`gmx mindist` is especially useful to find water
    molecules close to a region of interest.

:ref:`gmx do_dssp`
    Use the DSSP_ algorithm [Kabsch1983]_ to analyze the secondary structure
    (helices, sheets, ...). 


   
    
.. seealso::

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

     

