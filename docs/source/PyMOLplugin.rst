PyMOL plugin
============

.. _ConservedWaterSearch examples: https://conservedwatersearch.readthedocs.io/en/latest/examples.html
.. _ConservedWaterSearch documentation: https://conservedwatersearch.readthedocs.io/en/latest/conservedwaters.html

This plugin provides tools for identification, classification and analysis of
conserved water molecules and their networks from molecular dynamics
simulation trajectories.

The main interface consists of the following sections:

Input Options
-------------

Input trajectories **MUST** be aligned and centered before analysis. The
plugin will not perform any alignment or centering. If the system is
centered on the protein, it usually works well. However, in instances
where the area of interest is close to the edge of the simulation box,
system should be centered on the position of the center of mass of the area of
interest.

**Selection Input Type**: Choose the type of input you want to provide. Options include:

.. tabs::

  .. tab:: **PyMOL Selection**

    A PyMOL selection is first exported to a PDB file and water
    molecules are extracted via the water selection panel with
    MDAnalysis.

    **PyMOL Selection**: Enter the name of the PyMOL object or a
    selection string which contains waters (can contain other elements
    as well) in your PyMOL window. It is recommended that the selection
    includes both the system and the water molecules of interest.

    **Test Input**: Click this button to test the provided input.

  .. tab:: **From Files**

    Provide trajectory and topology files.

    **Trajectory and Topology Files**: Browse and select the appropriate
    files for your analysis. For list of supported file types check
    MDAnalysis documentation `here
    <https://www.mdanalysis.org/docs/documentation_pages/coordinates/init.html#supported-coordinate-formats>`_.
    File formats that can contain both trajectory and topology data such
    as PDB, should only be provided using trajectory field.

    **Every Nth Frame**: Specify how many frames to skip for analysis.

    **Test Input**: Click this button to test the provided input.

  .. tab:: **CWS Input**

    Use a pre-existing ConservedWaterSearch (CWS) input file. Provide the clustering data set and specify the number of frames.

    **CWS input data file**: Browse and select the CWS input data file.
    This data file should contain the clustering data set. If only
    oxygens are clustered every row should contain the x, y, and z
    coordinates of oxygen atoms for every frame to be studied. For
    analysis of water types (recommended) each row should contain the
    x, y, and z coordinates of oxygen and hydrogen atoms. For more
    information see `ConservedWaterSearch examples`_.
    A CWS input data file can be generated when providing a PyMOL
    selection or trajectory and topology files by using the **Export CWS
    Input Data** button in the **Water Selection** section.

    **Number of snapshots**: Number of snapshots in the trajectory used
    to generate the CWS input data file.

  .. tab:: **Restart**

    Restart a previous water clustering analysis. Provide the partial
    results file and CWS restart data file for restarting the analysis.

    **Partial Results File**: Browse and select the partial results file
    for restarting the analysis. This file should contain the partial
    results of the clustering analysis. For more information see
    `ConservedWaterSearch examples`_.

    **CWS restart data file**: Browse and select the CWS restart data
    file. This file has the same structure as CWS input data file.
    This data file should contain the clustering data set. If only
    oxygens are clustered every row should contain the x, y, and z
    coordinates of oxygen atoms for every frame to be studied. For
    analysis of water types (recommended) each row should contain the
    x, y, and z coordinates of oxygen and hydrogen atoms. For more
    information see `ConservedWaterSearch examples`_.

Water Selection
---------------
This section provides guidance on selecting water molecules for
analysis. Proper water selection ensures that the analysis focuses on
relevant waters in your structure. Users should take great care when
choosing the waters of interest. Making an appropriate choice is
crucial, as selecting a large number of water molecules can
significantly slow down the analysis process. Typically, the most
relevant water molecules are those situated near the binding site,
proteins, or other surfaces. While you can analyze bulk waters, they
tend to diffuse freely. If results are produced from these waters,
especially those classified as WCW type, exercise caution when
interpreting the data.

- **Water Selection Center**: Choose the method for selecting the center
  of water molecules selection. Options include:

.. tabs::

  .. tab:: **Geometric Mean**

    Provide a MDAnalysis selection string. The geometric mean of the
    selected atoms will be used as the center. This selection is handled
    by MDAnalysis whose selection language is similar in most instances
    but not identical to PyMOL. For more information see `MDanalysis
    <https://www.mdanalysis.org/docs/documentation_pages/selections.html#simple-selections>`_.

  .. tab:: **XYZ**

    Specify the x, y, and z coordinates for center of water selection.


.. tabs::
  .. tab:: **Key Residue and Atom Names**

    Users can provide residue name assigned to water molecules, as well
    as the atom names for oxygen and hydrogen atoms in water molecules.
    Alternatively, the plugin offers an automatic option, which attempts
    to identify water residue names and atom names using conventions
    from widely-used MD engines and tools.

    **Solvent Residue Name**: Specify the name of the solvent residue or
    opt for automatic detection.

    **Water Oxygen Atom Name**: Specify the name of the water oxygen or
    opt for automatic detection.

    **Water Hydrogen Atom Name**: Specify the name of the water hydrogen or
    opt for automatic detection.

**Distance**: Specify the distance from the center for water selection
inside which waters shall be selected for analysis.

.. tabs::
  .. tab:: Buttons

    **Test Selection**: Click this button to test the water selection.

    **Export CWS Input Data**: Click this button to export the CWS input data to a file.

Compute results
---------------

.. tabs::

   .. tab:: Water Clustering

     Compute conserved waters and classify them into Fully, Half or
     Weakly Conserved water molecules. More information can be found in
     the `ConservedWaterSearch documentation`_.

     - **Clustering Method**. Choose the clustering method. Options include:

     .. tabs::

        .. tab:: **QMSRC**

          Quick Multi-Stage Re-Clustering procedure.
          The best ratio of quality and speed.

        .. tab:: **MSRC**

          Multi-Stage Re-Clustering procedure.
          Very slow, but very accurate.

        .. tab:: **SC**

          Single Clustering.
          Very fast, but not very accurate. Might work well for buried
          binding sites.

     - **Clustering Algorithm**. Choose the clustering algorithm.
       Options include:

     .. tabs::

        .. tab:: **HDBSCAN**

          Faster, but produces slightly worse clusters.

        .. tab:: **OPTICS**

          Slightly slower, but produces slightly better clusters.

     **Water Types for Clustering**: Select the types of water molecules
     for clustering. In principle users should choose ``FCW``, ``HCW`` and
     ``WCW``. In some cases it might make sense to leave ``WCW`` out. This
     will also reduce the time for the analysis by about a third. For more
     information see `ConservedWaterSearch documentation`_. ``OnlyO``
     option clusters only oxygen atoms without considering hydrogen
     orientation data.

     **Clustering Options**: Depending on the chosen method, provide the
     necessary parameters. Best to leave as is. For large number of snapshots
     (>1000) it is recommended to increase the value of ``EveryMinsamp`` to
     not more than 10% of the number of snapshots (if using QMSRC or MSRC).

     **Compute Clustering Button**: Click this button to start the clustering
     analysis.

     **Advanced Settings**

     Users are discouraged to change the default values for the advanced
     settings, except for number of threads setting under ``njobs``.

     .. note:: Number of threads: using more thread than 1 will often decrease the performance of the clustering procedure. Only use more than 1 thread if you have a very large system (thousands of frames) or a large selection of water molecules per frame.

   .. tab:: Water Densty Map

      Computes oxygen density maps by binning the location of oxygen atoms to a
      3D grid. After the histogram has been computed, results are convolved with
      a gaussian whose width is defined by the oxygen van der Waals radius. Use
      the isomesh slider to adjust the isomesh value for the density map. The
      slider can also be used after the map has been computed.

      **Grid Bin (Delta)**: Specify the bin size for the density map.

      **Output File Name**: Specify the name of the output file for the density map.

      **Compute Density Map Button**: Click this button to calculate the oxygen density map.

      **Isomesh Value Slider**: Adjust the slider to change the isomesh value for the density map. It can be used after the density map was computed, and it will update the computed map.
