Water Network Analysis Plugin for PyMOL
=======================================

This plugin provides tools for analyzing conserved waters in biomolecular structures.

Main Interface
--------------

The main interface consists of the following sections:

Input Options
-------------

- **Selection Input Type**: Choose the type of input you want to provide. Options include:

  - **PyMOL Selection**: Directly select atoms or molecules in PyMOL.
  - **From Files**: Provide trajectory and topology files.
  - **CWS Input**: Use a pre-existing ConservedWaterSearch input.
  - **Restart**: Restart a previous analysis.

- **PyMOL Selection**: Enter the selection string for PyMOL.

  - **Test Input**: Click this button to test the provided input.

- **Trajectory and Topology Files**: Browse and select the appropriate files for your analysis.

  - **Every Nth Frame**: Specify which frames to analyze from the trajectory.

  - **Test Input**: Click this button to test the provided input.

- **Clustering Input**: Provide the clustering data set and specify the number of frames.

- **Restart Input**: Provide the partial results file and clustering data set for restarting the analysis.

Water Selection
---------------

- **Water Selection Center**: Choose the method for selecting the center of the water molecule. Options include:

  - **XYZ**: Specify the x, y, and z coordinates.
  - **Geometric Mean**: Provide a MDanalysis selection string. The geometric mean of the selected atoms will be used as the center.

- **Solvent Residue Name**: Specify the name of the solvent residue.

- **Water Oxygen Atom Name**: Specify the name of the water oxygen atom.

- **Water Hydrogen Atom Name**: Specify the name of the water hydrogen atom.

- **Distance**: Specify the distance from the center for water selection.

- **Test Selection**: Click this button to test the water selection.

- **Export CWS Input Data**: Click this button to export the CWS input data to a file.

Water Clustering Tab
--------------------

Compute conserved waters and classify them into several groups. More
information can be found in the `ConservedWaterSearch documentation <https://conservedwatersearch.readthedocs.io/en/latest/conservedwaters.html>`_.

Advanced Settings
-----------------

- Users are discouraged to change the default values for the advanced
  settings. Note on number of threads: often using more threads than 1
  will not improve the performance of the analysis, but slow it down.
  Only use more than 1 thread if you have a very large system (thousands
  of frames).

- **Clustering Method**: Choose the clustering method. Options include:

  - **QMSRC**: The best ratio of quality and speed.
  - **MSRC**: very slow, but very accurate.
  - **SC**: very fast, but not very accurate. Might work well for deep
    binding sites.

- **Clustering Algorithm**: Choose the clustering algorithm. Options include:

  - **OPTICS** - slightly slower, ubt produces slightly better clusters.
  - **HDBSCAN** - faster, but produces slightly worse clusters.

- **Water Types for Clustering**: Select the types of water molecules
  for clustering. In principle users should choose ``FCW``, ``HCW`` and
  ``WCW``. In some cases it might make sense to leave ``WCW`` out. This
  will also reduce the time for the analysis by about a third. For more
  information see `ConservedWaterSearch <https://conservedwatersearch.readthedocs.io/en/latest/conservedwaters.html>`_.

- **Clustering Options**: Depending on the chosen method, provide the
  necessary parameters. Best to leave as is. For large number of snapshots
  (>1000) it is recommended to increase the value of ``EveryMinsamp`` to
  not more than 10% of the number of snapshots. 

- **Compute Clustering**: Click this button to start the clustering analysis.

Water Density Map Tab
---------------------
Computes oxygen density maps by binning the location of oxygen atoms to
a 3D grid. Use the isomesh slider to adjust the isomesh value for the
density map. The slider can also be used after the map has been
computed.

- **Grid Bin (Delta)**: Specify the bin size for the density map.

- **Output File Name**: Specify the name of the output file for the density map.

- **Compute Density Map**: Click this button to calculate the oxygen density map.

- **Isomesh Value Slider**: Adjust the slider to change the isomesh value for the density map.

