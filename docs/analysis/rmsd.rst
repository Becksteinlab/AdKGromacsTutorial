.. -*- encoding: utf-8 -*-

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`



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

We compute the |Calpha| **RMSD** with `g_rms`_ with respect to the
reference starting structure (the one used for creating the :file:`md.tpr`
file). Work in a separate analysis directory::

  mkdir analysis/RMSD && cd analysis/RMSD

First we **create an index file** for the |Calpha| atoms
[#default_ndx_groups]_. Use `make_ndx`_ to create a file
:file:`ca.ndx` that contains the |Calpha| atoms as an *index
group*. Start :program:`make_ndx` and use :file:`md.tpr` as input; the
output index file will be :file:`ca.ndx`::

  gmx make_ndx -f ../../MD/md.tpr -o CA.ndx

Use `make_ndx`_ interactively by typing the following commands [#scripted_make_ndx]_::

  keep 1
  a CA
  name 1 Calpha
  q

(This sequence of commands only retains the "Protein" default
selection, then selects all atoms named "CA", renames the newly
created group to "Calpha", and saves and exits.)

You can look at :file:`CA.ndx` and see all the index numbers listed
under the heading ``[ Calpha ]``.

Run :program:`g_rms`, using our newly defined group as the selection
to fit and to compute the RMSD::

  printf "Calpha\nCalpha\n" | gmx rms -s ../../MD/md.tpr -f ../../MD/md.xtc -n CA.ndx -o rmsd.xvg -fit rot+trans

Note that the units are nm.

Plot the ``rmsd.xvg`` file as usual (you might want to use :kbd:`gmx rms
-xvg none` if you are processing with NumPy/matplotlib_):

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy
   t,rmsd = numpy.loadtxt("rmsd.xvg", unpack=True)
   fig = plt.figure(figsize=(5,2.5)) 
   ax = fig.add_subplot(111)
   fig.subplots_adjust(bottom=0.2)
   ax.set_xlabel("time $t$ (ps)")
   ax.set_ylabel(r"C$_\alpha$ RMSD (nm)")
   ax.fill_between(t,rmsd, color="blue", linestyle="-", alpha=0.1) 
   ax.plot(t,rmsd, color="blue", linestyle="-") 
   fig.savefig("rmsd_ca.png", dpi=300)
   fig.savefig("rmsd_ca.svg")
   fig.savefig("rmsd_ca.pdf")

.. figure:: /figs/rmsd_ca.*
   :scale: 80%
   :alt: RMSD timeseries
   
   Root mean square distance (RMSD) of the |Calpha| atoms of AdK from
   the initial simulation frame.
   

   
.. links
.. _`g_rms`: http://manual.gromacs.org/current/online/g_rms.html
.. _`make_ndx`: http://manual.gromacs.org/current/online/make_ndx.html
.. _matplotlib: https://matplotlib.org

.. rubric:: Footnotes
	    
.. [#default_ndx_groups] Actually, we don't need to create the index
   group for the |Calpha| atoms ourselves because Gromacs
   automatically creates the group "C-alpha" as one of many default
   groups (other are "Protein", "Protein-H" (only protein heavy
   atoms), "Backbone" (N CA C), "Water", "non-Protein" (i.e. water and
   ions in our case but could also contain other groups such as drug
   molecule or a lipid membrane in more complicated simulations),
   "Water_and_ions". You can see these index groups if you just run
   :program:`make_ndx` on an input structure or if you interactively
   select groups in :program:`trjconv`, :program:`g_rms`, ...

   However, making the "Calpha" group yourself is a good exercise
   because in many cases there are no default index groups for the
   analysis you might want to do.

.. [#scripted_make_ndx] In scripts you can pipe all the interactive
   commands to `make_ndx`_ by using the :code:`printf ... | make_ndx`
   trick::
     
     printf "keep 0\ndel 0\na CA\nname 0 Calpha\nq\n" | make_ndx -f ../../MD/md.tpr -o CA.ndx

   This will accomplish the same thing as the interactive use
   described above.

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

