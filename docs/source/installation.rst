Installation
============
The easiest ways to install **WaterNetworkAnalysis** is using :code:`conda` from conda-forge:

.. code:: bash

    conda install -c conda-forge WaterNetworkAnalysis

Alternatively, WNA is also available on PyPi via :code:`pip`. However, because WNA depends on `ConservedWaterSearch <https://conservedwatersearch.readthedocs.io/en/latest/installation.html>`_ which requires hdbscan whose PyPi installation requires a C++ compiler (see `here <https://conservedwatersearch.readthedocs.io/en/latest/installation.html>`_ for more information) aditional dependencies have to be installed:

.. code:: bash

   conda install -c conda-forge cxx-compiler

`Pymol <https://pymol.org/2/>`_ is an optional dependency for visualisation and is not present on PyPi, however WNA can be installed and used without it (bar pymol visualisation features). Pymol can be installed using :code:`conda`:

.. code:: bash

   conda install -c conda-forge pymol-open-source
 
Finally, to install WNA via :code:`pip` use:

.. code:: bash

   pip install WaterNetworkAnalysis

For more information on CWS dependencies also see `CWS installation guide <https://conservedwatersearch.readthedocs.io/en/latest/installation.html>`_.

PyMOL plugin installation
=========================

Known Issues with dependencies
==============================

:code:`AttributeError: 'super' object has no attribute '_ipython_display_'`
Some versions of Jupyter notebook are incpompatible with ipython (`see here <https://stackoverflow.com/questions/74279848/nglview-installed-but-will-not-import-inside-juypter-notebook-via-anaconda-navig>`_). To resolve install version of :code:`ipywidgets<8` using :code:`conda`: 

.. code:: bash

   conda install "ipywidgets <8" -c conda-forge

or :code:`pip`:

.. code:: bash

   pip install ipywidgets==7.6.0