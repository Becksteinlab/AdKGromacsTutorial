.. -*- encoding: utf-8 -*-

.. include:: /includes/defs.rst
.. include:: /includes/links.rst


.. _distances:

===========
 Distances
===========

Distances can be a very useful observable to track conformational
changes. They can often be directly related to real experimental
observables such as NOEs from NMR experiments or distances from
cross-linking or FRET experiments.

Here we calculate a simple distance between two |Calpha| atoms as an
approximation to the distance between two chromophores attached to the
correspoding residues isoleucine 52 (*I52*) and lysine 145 (*K145*) used
in a FRET experiment [Henzler-Wildman2007]_.

First we need to create an index file containing the two groups::

 mkdir -p analysis/dist/I52_K145 && cd analysis/dist/I52_K145
 gmx make_ndx -f ../../../MD/md.tpr -o I52_K145.ndx

Use interactive commands like the following [#ndx_selections]_ ::

 keep 0
 del 0
 r 52 & a CA
 name 0 I52
 r 145 & a CA
 name 1 K145
 q

to generate the index file :file:`I52_K145.ndx`.
 
Then run :ref:`gmx distance` and compute the distance between the two atoms::

  printf "I52\nK145\n" | gmx distance -s ../../../MD/md.tpr -f ../../../MD/md.xtc -n I52_K145.ndx -o dist.xvg

The :file:`dist.xvg` file contains the distance in nm for each time
step in ps, which can be plotted [#plot_distance]_.

.. figure:: /figs/d_I52_K145_ca.*
   :scale: 80%
   :alt: Timeseries of the distance between |Calpha| atoms of I52 and K145.
   
   Timeseries of the distance between the |Calpha| atoms of I52 (NMP
   domain) and K145 (LID domain).

(You can also use the centered and fitted trajectory
:file:`md_fit.xtc` as an input instead of :file:`md.xtc` to make sure
that the distance calculation does not contain any jumps due to
periodic boundary effects, or use :ref:`gmx mindist`.)

.. SeeAlso:: [Beckstein2009]_ for a discussion of FRET distances in AdK.



.. rubric:: Footnotes
	    
.. [#ndx_selections] Note that one has to be careful when selecting
   residue ids in :program:`make_ndx`. It is often the case that a PDB
   file does not contain all residues, e.g. residues 1--8 might be
   unresolved in the experiment and thus are missing from the PDB
   file. The file then simply starts with residue number 9. Gromacs,
   however, typically *renumbers residues so that they start at
   1*. Thus, in this hypothetical case, a residue that might be
   referred to in the literature as "residue 100" might actually be
   residue 92 in the simulation (:math:`N^\mathrm{sim}_\mathrm{res} =
   N^\mathrm{PDB}_\mathrm{res} - (\mathrm{min}
   N^\mathrm{PDB}_\mathrm{res} - 1)`). Thus, if you wanted to select
   the |Calpha| atom of residue 100 you would need to select :kbd:`r
   92 & a CA` in :program:`make_ndx`.
	     
.. [#plot_distance]

   To plot in Python, make sure to *not* write xvg legend information
   to the output file (using the ``-xvg none`` flag)

   .. code-block:: bash

      printf "I52\nK145\n" | \
           gmx distance -s ../../../MD/md.tpr -f ../../../MD/md.xtc \
	                -n I52_K145.ndx -o dist.xvg -xvg none
		   
   so that you can easily read the data with :func:`numpy.loadtxt`:

   .. code-block:: python
   
      import matplotlib.pyplot as plt
      import numpy
      
      t,d,x,y,z = numpy.loadtxt("dist.xvg", unpack=True)
      
      fig = plt.figure(figsize=(5,2.5)) 
      ax = fig.add_subplot(111)
      fig.subplots_adjust(bottom=0.2)
      
      ax.fill_between(t, d, color="orange", linestyle="-", alpha=0.1) 
      ax.plot(t, d, color="orange", linestyle="-", label="I52-K145")

      ax.set_xlabel("time $t$ (ps)")
      ax.set_ylabel(r"C$_\alpha$ distance (nm)")
      ax.legend(loc="best")
      
      fig.savefig("d_I52_K145_ca.png", dpi=300)
      fig.savefig("d_I52_K145_ca.svg")
      fig.savefig("d_I52_K145_ca.pdf")

