.. -*- encoding: utf-8 -*-

.. |kJ/mol/nm**2| replace:: kJ mol\ :sup:`-1` nm\ :sup:`-2`
.. |Calpha| replace:: C\ :sub:`α`

.. _RMSF:

======
 RMSF
======

The residue root mean square fluctuation **RMSF** is a measure of the
flexibility of a residue. It is typically calculated for the |Calpha|
atom of each residue and is then simply the square root of the
variance of the fluctuation around the average position:

.. math::

   \rho^{\mathrm{RMSF}}_i = \sqrt{\left\langle 
        \left(\mathbf{r}_i - \langle \mathbf{r}_i \rangle \right)^2 
        \right\rangle}

Use the :file:`CA.ndx` file from the :ref:`RMSD` calculation with
`gmx rmsf`_ ::

  mkdir analysis/RMSF && cd analysis/RMSF
  printf "Calpha\n" | g_rmsf -s ../../MD/md.tpr -f ../../MD/md.xtc -n ../RMSD/CA.ndx -o rmsf.xvg -fit

A plot [#plot_rmsf]_ of :math:`\rho^{\mathrm{RMSF}}_{i}` versus residue number *i*
shows the regions of high flexibility as peaks in the plot. Note that
a 100-ps simulation might be too short to obtain a meaningful RMSF
profile. 

.. figure:: /figs/rmsf_ca.*
   :scale: 80%
   :alt: Per-residue RMSF
   
   Root mean square fluctuation (RMSF) of the |Calpha| atoms of AdK.





.. rubric:: Comparison with crystallographic B-factors

You can compare the RMSF to the isotropic atomic crystallographic
B-factors, which are related by [Willis1975]_

.. math::

   B_{i} = \frac{8\pi^2}{3} (\rho^{\mathrm{RMSF}}_{i})^2

(In this case you would want to calculate the RMSF for all heavy
(i.e. non-hydrogen) atoms. You don't need to build and use a separate
index file file: simply choose the default group "Protein-H" ("protein
without hydrogens")).

.. Note:: Gromacs RMSF are in units of nm and B-factors are
          typically measured in Å\ :sup:`2`.

It is straightforward to write Python code that calculates the
B-factor from the RMSF in :file:`rmsf.xvg` and it is also easy to
extract the B-factor ("temperatureFactor") from columns 61-66 in the
`ATOM record of a PDB file`_.


.. _`gmx rmsf`: http://manual.gromacs.org/documentation/current/onlinehelp/gmx-rmsf.html
.. _`ATOM record of a PDB file`:
   http://www.wwpdb.org/documentation/format33/sect9.html#ATOM

.. rubric:: Footnotes

.. [#plot_rmsf]

   To plot in Python, make sure to *not* write xvg legend information
   to the output file (using the ``-xvg none`` flag)

   .. code-block:: bash

      printf "Calpha\n" | \
          g_rmsf -s ../../MD/md.tpr -f ../../MD/md.xtc -n  ../RMSD/CA.ndx \
                 -o rmsf.xvg -fit -xvg none
		   

   so that you can easily read the data with :func:`numpy.loadtxt`:

   .. code-block:: python
		   
      import matplotlib.pyplot as plt
      import numpy as np
      
      resid, rmsf = np.loadtxt("rmsf.xvg", unpack=True)
      
      fig = plt.figure(figsize=(5,2.5))
      ax = fig.add_subplot(111)
      fig.subplots_adjust(bottom=0.2)
      ax.set_xlabel("residue number")
      ax.set_ylabel(r"C$_\alpha$ RMSF (nm)")
      ax.fill_between(resid, rmsf, color="red", linestyle="-", alpha=0.1)
      ax.plot(resid, rmsf, color="red", linestyle="-")
      ax.set_xlim(resid.min(), resid.max())
      fig.savefig("rmsf_ca.png", dpi=300)
      fig.savefig("rmsf_ca.svg")
      fig.savefig("rmsf_ca.pdf")
	    
