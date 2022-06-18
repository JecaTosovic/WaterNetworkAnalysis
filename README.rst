WaterNetworkAnalysis
==============================
.. image:: https://readthedocs.org/projects/waternetworkanalysis/badge/?version=latest
    :target: https://waternetworkanalysis.readthedocs.io/en/latest/?badge=latest
.. image:: https://badge.fury.io/py/WaterNetworkAnalysis.svg
    :target: https://badge.fury.io/py/WaterNetworkAnalysis

The WaterNetworkAnalysis (WNA) Python package serves as a set of tools for input preparation for `ConservedWaterSearch <https://github.com/JecaTosovic/ConservedWaterSearch>`_ python package from MD trajectories and their analysis.
Using the ConservedWaterSearch package the conserved water molecules can be further classified into 3 distinct conserved water types based on their hydrogen orientation: Fully Conserved Waters (FCW), Half Conserved Waters (HCW) and Weakly conserved waters (WCW) - see figure below for examples.
WNA can be used to create PyMol or nglview visualisations of conserved water networks for drug discovery or materials science purposes.

.. image:: figs/WaterTypes.png
  :width: 600

Important links
=================
	- `Documentation <https://waternetworkanalysis.readthedocs.io/en/latest/>`_: hosted on Read The Docs
	- `GitHub repository <https://github.com/JecaTosovic/WaterNetworkAnalysis>`_: source code/contribute code
	- `Issue tracker <https://github.com/JecaTosovic/WaterNetworkAnalysis/issues>`_: Report issues/ request features

Related Tools
=================
	- `ConservedWaterSearch <https://github.com/JecaTosovic/ConservedWaterSearch>`_: prepare trajectories  and analyse results for/from conservedwatersearch

Citation
===============
Coming soon.

Installation
===============
The easiest ways to install **WaterNetworkAnalysis** is using pip:

.. code:: bash

   pip install WaterNetworkAnalysis

Conda builds will be available soon.


Example
===============
The following example shows how to use **WaterNetworkAnalysis** to prepare a MD trajectory and analyse the results for determination of conserved water networks.

.. code:: python

   import WaterNetworkAnalysis as WNA
   # load some example
   # Run classification
   # Create visualisations
   # plot results


