Water Network Analysis Plugin for PyMOL
=======================================

This plugin provides tools for analyzing water networks in molecular structures.

Main Interface
--------------

The main interface consists of the following sections:

- **Top Buttons**:
  
  - **Instructions**: Click this button to open the help window.
  - **Bibliography**: Click this button to view the citations related to this plugin.

Input Options
-------------

- **Selection Input Type**: Choose the type of input you want to provide. Options include:

  - **PyMOL Selection**: Directly select atoms or molecules in PyMOL.
  - **From Files**: Provide trajectory and topology files.
  - **CWS Input**: Use a pre-existing CWS input.
  - **Restart**: Restart a previous analysis.

- **PyMOL Selection**: Enter the selection string for PyMOL.

- **Trajectory and Topology Files**: Browse and select the appropriate files for your analysis.

- **Every Nth Frame**: Specify which frames to analyze from the trajectory.

- **Test Input**: Click this button to test the provided input.

- **Clustering Input**: Provide the clustering data set and specify the number of frames.

- **Restart Input**: Provide the partial results file and clustering data set for restarting the analysis.

Water Selection
---------------

- **Water Selection Center**: Choose the method for selecting the center of the water molecule. Options include:

  - **XYZ**: Specify the x, y, and z coordinates.
  - **Geometric Mean**: Provide a PyMOL selection string.

- **Solvent Residue Name**: Specify the name of the solvent residue.

- **Water Oxygen Atom Name**: Specify the name of the water oxygen atom.

- **Water Hydrogen Atom Name**: Specify the name of the water hydrogen atom.

- **Distance**: Specify the distance for water selection.

- **Test Selection**: Click this button to test the water selection.

- **Export CWS Input Data**: Click this button to export the CWS input data to a file.

Advanced Settings (if allowed)
------------------------------

- **Aligning Method**: Choose the method for aligning the structure.

Water Clustering Tab
--------------------

- **Clustering Method**: Choose the clustering method. Options include:

  - **QMSRC**
  - **MSRC**
  - **SC**

- **Clustering Algorithm**: Choose the clustering algorithm. Options include:

  - **OPTICS**
  - **HDBSCAN**

- **Water Types for Clustering**: Select the types of water molecules for clustering.

- **Clustering Options**: Depending on the chosen method, provide the necessary parameters.

- **Compute Clustering**: Click this button to start the clustering analysis.

Water Density Map Tab
---------------------

- **Grid Bin (Delta)**: Specify the bin size for the density map.

- **Output File Name**: Specify the name of the output file for the density map.

- **Compute Density Map**: Click this button to calculate the oxygen density map.

- **Isomesh Value Slider**: Adjust the slider to change the isomesh value for the density map.

