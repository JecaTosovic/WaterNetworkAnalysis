Installation
============

WaterNetworkAnalysis can be installed as a python package or as a PyMOL
plugin.

.. include:: ../../README.rst
    :start-line: 29
    :end-before: Installation of the PyMOL plugin

Python package installation
---------------------------

The easiest ways to install **WaterNetworkAnalysis** is using :code:`conda` from conda-forge:

.. code:: bash

    conda install -c conda-forge WaterNetworkAnalysis

Alternatively, WNA is also available on PyPi via :code:`pip`:

.. code:: bash

   pip install WaterNetworkAnalysis

`Pymol <https://pymol.org/2/>`_ is an optional dependency for visualisation and is not present on PyPi, however WNA can be installed and used without it (bar pymol visualisation features). Pymol can be installed using :code:`conda`:

.. code:: bash

   conda install -c conda-forge pymol-open-source
 
For more information on CWS dependencies also see `CWS installation guide <https://conservedwatersearch.readthedocs.io/en/latest/installation.html>`_.

PyMOL plugin installation
-------------------------

This guide provides detailed installation instructions for Linux, Mac,
and Windows users. We recommend using `conda` or `mamba` to create a new
environment with at least Python 3.9.

Note that PyMOL which can be downloaded from the `PyMOL website
<https://pymol.org/2/>`_ comes with python 3.7 which is not supported by 
WaterNetworkAnalysis or ConservedWaterSearch. For this reason users
will have to install `mamba/conda` and create a new environment with
python version greater or equal to 3.9 and install PyMOL in that environment.
The plugin has been tested with PyMOL version > 2.5.0.

Prerequisites
.............

- Ensure you have ``conda`` or ``mamba`` installed. If not, download and
  install `miniforge <https://conda-forge.org/miniforge/>`_ , or
  `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or
  `Anaconda <https://www.anaconda.com/products/distribution>`_.

- For users who wish to use the paid version of PyMOL, ensure you have a
  valid license.

- Once ``conda/mamba`` has been installed, make sure you activate the
  installation by either sourcing your ``.bashrc`` or restarting your terminal.

Installation Steps
..................

1. **Create a new conda environment (use mamba instead of conda if you
   opted for mamba)**:

     .. code-block:: bash

        conda create -n myenv python=3.9

   Replace ``myenv`` with your preferred environment name.

2. **Activate the environment**:

   - **Linux & Mac**:

     .. code-block:: bash

        conda activate myenv

   - **Windows**:

     .. code-block:: bash

        activate myenv

3. **Install PyMOL**:

   - **Open-source version**:

     .. code-block:: bash

        conda install -c conda-forge pymol-open-source


   - **Paid version**:

     .. code-block:: bash

        conda install -c schrodinger pymol-bundle


   macOS users may need to install the extra packages. For more
   information see `PyMOL documentation <https://pymol.org/conda/>`_. To
   test if the installation was successful users should just be able to
   type the following in their terminal:
   
     .. code-block:: bash
    
        pymol

   
   Users with a license should download their license file from the
   `PyMOL website <https://pymol.org/2/>`_ and activate it by going to
   Help -> Install new License File in main PyMOL window.

4. **Install dependencies**:

   WaterNetworkAnalysis is the main dependency and can be installed via:

     .. code-block:: bash

        conda install -c conda-forge WaterNetworkAnalysis


6. **Install the WaterNetworkAnalysis plugin in PyMOL**:

   The plugin is a single file located `here
   <https://github.com/JecaTosovic/WNA_PyMOL_plugin>`_. In PyMOL, go to
   `Plugin > Plugin Manager > Install New Plugin` and select the
   ``WNA_PyMOL_plugin.py`` file from the ``WNA_PyMOL`` folder. The
   plugin can then be accessed from the plugin drop-down menu.

Troubleshooting
---------------

- If you encounter any issues, ensure you're using the correct Python
  version and that all packages are installed with their specified
  versions.

- For further assistance, refer to the official documentation or contact
  the support team.


Known Issues with dependencies
------------------------------

:code:`AttributeError: 'super' object has no attribute '_ipython_display_'`
Some versions of Jupyter notebook are incpompatible with ipython (`see here <https://stackoverflow.com/questions/74279848/nglview-installed-but-will-not-import-inside-juypter-notebook-via-anaconda-navig>`_). To resolve install version of :code:`ipywidgets<8` using :code:`conda`: 

.. code:: bash

   conda install "ipywidgets <8" -c conda-forge

or :code:`pip`:

.. code:: bash

   pip install ipywidgets==7.6.0