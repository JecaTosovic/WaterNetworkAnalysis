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

This guide provides detailed installation instructions for Linux, Mac,
and Windows users. We recommend using `conda` or `mamba` to create a new
environment with at least Python 3.9.

Note that PyMOL which can be downloaded from the `PyMOL website
<https://pymol.org/2/>`_ comes with python 3.7 which is not supported by
WaterNetworkAnalysis and ConservedWaterSearch. For this reason users
will have to install mamba/conda and create a new environment with
python version greater or equal to 3.9 and install PyMOL in that environment.
The plugin has been tested with PyMOL 2.5.0.

Prerequisites
-------------

- Ensure you have `conda` or `mamba` installed. If not, download and
  install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or
  [Anaconda](https://www.anaconda.com/products/distribution).

- For users who wish to use the paid version of PyMOL, ensure you have a
  valid license. Follow the instructions provided by the PyMOL vendor to
  activate your license.

- Once conda/mamba has been installed, make sure you activate the
  installation by either sourcing your .bashrc or restarting your terminal.

Installation Steps
------------------

1. **Create a new conda environment**:

   .. code-block:: bash

      conda create -n myenv python=3.9

   Replace `myenv` with your preferred environment name.

2. **Activate the environment**:

   - **Linux & Mac**:

     .. code-block:: bash

        conda activate myenv

   - **Windows**:

     .. code-block:: bash

        activate myenv

3. **Install the required packages**:

   Use the following command to install the necessary packages:

   .. code-block:: bash

      conda install -c conda-forge pymol-open-source
      ConservedWaterSearch WaterNetworkAnalysis


4. **Install PyMOL**:

   - **Open-source version**:

     .. code-block:: bash

        conda install -c conda-forge pymol

   - **Paid version**:

     Follow the installation instructions provided by the PyMOL vendor.

5. **Verify Installation**:

   After installing all the required packages, you can verify the PyMOL
   installation by simply starting it from the terminal:

   .. code-block:: bash

      pymol


6. **Install the WaterNetworkAnalysis plugin in PyMOL**:

   In PyMOL, go to `Plugin > Plugin Manager > Install New Plugin` and
    select the `WaterNetworkAnalysis.py` file from the
    `WaterNetworkAnalysis` folder in the installation directory.

Troubleshooting
---------------

- If you encounter any issues, ensure you're using the correct Python
  version and that all packages are installed with their specified
  versions.

- For further assistance, refer to the official documentation or contact
  the support team.


Known Issues with dependencies
==============================

:code:`AttributeError: 'super' object has no attribute '_ipython_display_'`
Some versions of Jupyter notebook are incpompatible with ipython (`see here <https://stackoverflow.com/questions/74279848/nglview-installed-but-will-not-import-inside-juypter-notebook-via-anaconda-navig>`_). To resolve install version of :code:`ipywidgets<8` using :code:`conda`: 

.. code:: bash

   conda install "ipywidgets <8" -c conda-forge

or :code:`pip`:

.. code:: bash

   pip install ipywidgets==7.6.0