.. -*- encoding: utf-8 -*-

.. include:: /includes/defs.rst
.. include:: /includes/links.rst

	     
.. _radius-of-gyration:

====================
 Radius of gyration
====================

The radius of gyration measures the compactness of a protein
structure.

.. math::

   R_\mathrm{gyr}^2 = \frac{1}{M}\sum_{i=1}^{N} m_i(\mathbf{r}_i - \mathbf{R})^2

where :math:`M = \sum_{i=1}^{N} m_i` is the total mass and
:math:`\mathbf{R} = N^{-1}\sum_{i=1}^{N} \mathbf{r}_i` is the center of
mass of the protein consisting of :math:`N` atoms.

The Gromacs tool :ref:`gmx gyrate` can be used to compute the radius of
gyration for the whole protein (using the pre-defined "Protein" index
group) ::

  mkdir analysis/rgyr && cd analysis/rgyr
  echo Protein | gmx gyrate -s ../../MD/md.tpr -f ../../MD/md.xtc -o gyrate.xvg
  
and the resulting time series in file :file:`gyrate.xvg` can be
plotted [#plot_rgyr]_.

.. figure:: /figs/rgyr.*
   :scale: 80%
   :alt: Timeseries of the radius of gyration.
   
   Timeseries of the radius of gyration computed for the whole
   protein.

(In principle, :math:`R_\mathrm{gyr}` can indicate changes in
conformation but the simulation time in our test trajectory is too
short to reveal the large conformational transition that AdK can
undergo [Beckstein2009]_.)
	 


.. rubric:: Footnotes

.. [#plot_rgyr]

   To plot in Python, make sure to *not* write xvg legend information
   to the output file (using the ``-xvg none`` flag)

   .. code-block:: bash

      echo Protein | \
          gmx gyrate -s ../../MD/md.tpr -f ../../MD/md.xtc -o gyrate.xvg -xvg none

   so that you can easily read the data with :func:`numpy.loadtxt`:

   .. code-block:: python   

      import matplotlib.pyplot as plt
      import numpy
      
      t,data,x,y,z = numpy.loadtxt("gyrate.xvg", unpack=True)
      
      fig = plt.figure(figsize=(5,2.5))
      ax = fig.add_subplot(111)
      fig.subplots_adjust(bottom=0.2)
      
      ax.fill_between(t,data, color="magenta", linestyle="-", alpha=0.1)
      ax.plot(t,data, color="magenta", linestyle="-")

      ax.set_xlabel("time $t$ (ps)")
      ax.set_ylabel(r"protein $R_\mathrm{gyr}$ (nm)")
      
      fig.savefig("rgyr.png", dpi=300)
      fig.savefig("rgyr.svg")
      fig.savefig("rgyr.pdf")
