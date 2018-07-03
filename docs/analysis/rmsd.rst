.. -*- encoding: utf-8 -*-

.. include:: /includes/defs.rst
.. include:: /includes/links.rst


.. _RMSD:

======
 RMSD
======

The RMSD is the root mean squared Euclidean distance in *3N*
configuration space as function of the time step,

.. math::

   \rho^{\mathrm{RMSD}}(t) = \sqrt{\frac{1}{N} \sum_{i=1}^{N}\left(\mathbf{r}_{i}(t) - \mathbf{r}_{i}^{\mathrm{ref}}\right)^2}

between the current coordinates :math:`\mathbf{r}_{i}(t)` at time *t*
and the reference coordinates :math:`\mathbf{r}_{i}^{\mathrm{ref}}`.

We compute the |Calpha| **RMSD** with :ref:`gmx rms` with respect to the
reference starting structure (the one used for creating the :file:`md.tpr`
file). Work in a separate analysis directory::

  mkdir analysis/RMSD && cd analysis/RMSD

First we **create an index file** for the |Calpha| atoms
[#default_ndx_groups]_. Use :ref:`gmx make_ndx` to create a file
:file:`ca.ndx` that contains the |Calpha| atoms as an *index
group*. Start :program:`make_ndx` and use :file:`md.tpr` as input; the
output index file will be :file:`ca.ndx`::

  gmx make_ndx -f ../../MD/md.tpr -o CA.ndx

Use :ref:`gmx make_ndx` interactively by typing the following commands
[#scripted_make_ndx]_::

  keep 1
  a CA
  name 1 Calpha
  q

(This sequence of commands only retains the "Protein" default
selection, then selects all atoms named "CA", renames the newly
created group to "Calpha", and saves and exits.)

You can look at :file:`CA.ndx` and see all the index numbers listed
under the heading ``[ Calpha ]``.

Run :program:`gmx rms`, using our newly defined group as the selection
to fit and to compute the RMSD::

  printf "Calpha\nCalpha\n" | gmx rms -s ../../MD/md.tpr -f ../../MD/md.xtc -n CA.ndx -o rmsd.xvg -fit rot+trans

Note that the units are nm.

Plot the timeseries data in the :file:`rmsd.xvg` [#plot_rmsd]_.

.. figure:: /figs/rmsd_ca.*
   :scale: 80%
   :alt: RMSD timeseries
   
   Root mean square distance (RMSD) of the |Calpha| atoms of AdK from
   the initial simulation frame.
   

   
.. rubric:: Footnotes
	    
.. [#default_ndx_groups] Actually, we don't need to create the index
   group for the |Calpha| atoms ourselves because Gromacs
   automatically creates the group "C-alpha" as one of many default
   groups (other are "Protein", "Protein-H" (only protein heavy
   atoms), "Backbone" (N CA C), "Water", "non-Protein" (i.e. water and
   ions in our case but could also contain other groups such as drug
   molecule or a lipid membrane in more complicated simulations),
   "Water_and_ions". You can see these index groups if you just run
   :ref:`gmx make_ndx` on an input structure or if you interactively
   select groups in :ref:`gmx trjconv`, :ref:`gmx rms`, ...

   However, making the "Calpha" group yourself is a good exercise
   because in many cases there are no default index groups for the
   analysis you might want to do.

.. [#scripted_make_ndx] In scripts you can pipe all the interactive
   commands to :ref:`gmx make_ndx` by using the :code:`printf ... | gmx
   make_ndx` trick::
     
     printf "keep 0\ndel 0\na CA\nname 0 Calpha\nq\n" | \
          gmx make_ndx -f ../../MD/md.tpr -o CA.ndx

   This will accomplish the same thing as the interactive use
   described above.

.. [#plot_rmsd]

   If you use Python (namely NumPy_ and matplotlib_) then you might
   want to use :kbd:`gmx rms -xvg none` so that *no XVG legend
   information* is written to the output file

   .. code-block:: bash

      printf "Calpha\nCalpha\n" | \
          gmx rms -s ../../MD/md.tpr -f ../../MD/md.xtc -n CA.ndx \
	          -o rmsd.xvg -fit rot+trans -xvg none

   and you can easily read the data with :func:`numpy.loadtxt`:

   .. code-block:: python

      import matplotlib.pyplot as plt
      import numpy

      t,rmsd = numpy.loadtxt("rmsd.xvg", unpack=True)

      fig = plt.figure(figsize=(5,2.5)) 
      ax = fig.add_subplot(111)
      fig.subplots_adjust(bottom=0.2)

      ax.fill_between(t,rmsd, color="blue", linestyle="-", alpha=0.1) 
      ax.plot(t,rmsd, color="blue", linestyle="-")

      ax.set_xlabel("time $t$ (ps)")
      ax.set_ylabel(r"C$_\alpha$ RMSD (nm)")

      fig.savefig("rmsd_ca.png", dpi=300)
      fig.savefig("rmsd_ca.svg")
      fig.savefig("rmsd_ca.pdf")   
