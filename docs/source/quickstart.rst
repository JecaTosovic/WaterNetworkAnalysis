WaterNetworkAnalysis
====================
.. image:: https://readthedocs.org/projects/waternetworkanalysis/badge/?version=latest
    :target: https://waternetworkanalysis.readthedocs.io/en/latest/?badge=latest
.. image:: https://badge.fury.io/py/WaterNetworkAnalysis.svg
    :target: https://badge.fury.io/py/WaterNetworkAnalysis

The WaterNetworkAnalysis (WNA) Python package serves as a set of tools for input preparation for `ConservedWaterSearch <https://github.com/JecaTosovic/ConservedWaterSearch>`_ python package which identifies conserved water molecules from Molecular Dynamics (MD) trajectories.

.. image:: figs/Scheme.png
  :width: 600

WNA can also be used to create PyMol or nglview visualisations of conserved water networks for drug discovery or materials science purposes.

Important links
---------------
	- `Documentation <https://waternetworkanalysis.readthedocs.io/en/latest/>`_: hosted on Read The Docs
	- `GitHub repository <https://github.com/JecaTosovic/WaterNetworkAnalysis>`_: source code/contribute code
	- `Issue tracker <https://github.com/JecaTosovic/WaterNetworkAnalysis/issues>`_: Report issues/ request features

Related Tools
-------------
	- `ConservedWaterSearch <https://github.com/JecaTosovic/ConservedWaterSearch>`_: Analysis of conserved waters from simulation trajectories. For docs see `here <https://conservedwatersearch.readthedocs.io/en/latest/>`_.

Citation
--------
For citations and more infromation see `ConservedWaterSearch citation <https://conservedwatersearch.readthedocs.io/en/latest/citing.html>`_.

Installation
------------
The easiest ways to install **WaterNetworkAnalysis** is using pip:

.. code:: bash

   pip install WaterNetworkAnalysis

Conda builds will be available soon.

WNA depends on ConservedWaterSearch (CWS). To utilize full visualisation capabilities of CWS, `pymol <https://pymol.org/2/>`_ should be installed. Easiest way to install pymol is via ``conda``:

.. code:: bash

   conda install -c conda-forge pymol-open-source

For more information see `CWS installation guide <https://conservedwatersearch.readthedocs.io/en/latest/installation.html>`_.

Example
-------
The following example shows how to use **WaterNetworkAnalysis** to prepare a MD trajectory and analyse the results for determination of conserved water networks.

.. code:: python

   from WaterNetworkAnalysis import align_trajectory, get_center_of_selection,get_selection_string_from_resnums,extract_waters_from_trajectory
   from ConservedWaterSearch.water_clustering import WaterClustering
   from ConservedWaterSearch.utils import get_orientations_from_positions
   
   # MD trajectory filename
   trajectory="md.xtc"
   # topology filename
   topology="md.gro"
   # aligned trajectory filename
   alignedtrj = "aligned_trj.xtc"
   # aligned snapshot filename
   aligned_snap = "aligned.pdb"
   # distance to select water molecules around
   distance = 12.0
   # align the trajectory and save the alignment reference configuration
   align_trajectory(
       trajectory=trajectory,
       topology=topology,
       align_target_file_name=aligned_snap,
       output_trj_file=alignedtrj,
   )
   # define active site by aminoacid residue numbers
   active_site_resnums = [111, 112, 113, 122, 133, 138, 139, 142, 143, 157, 166, 167, 169, 170, 203, 231, 232, 238]
   # find centre of the active site in aligned trajectory
   selection_centre = get_center_of_selection(
       get_selection_string_from_resnums(active_site_resnums),
       trajectory=alignedtrj,
       topology=topology,
   )
   # extract water coordinates of interest around selection centre
   coordO, coordH =  extract_waters_from_trajectory(
       trajectory=alignedtrj,
       topology=topology,
       selection_center=selection_centre, 
       dist=distance
   )
   # start the clustering procedure
   Nsnaps = 200
   WC=WaterClustering(nsnaps= Nsnaps)
   # perform multi stage reclustering
   WC.multi_stage_reclustering(*get_orientations_from_positions(coordO,coordH))
   # visualise results with pymol
   WC.visualise_pymol(aligned_snap, active_site_ids=active_site_resnums, dist=distance)



.. image:: figs/Results.png
  :width: 600