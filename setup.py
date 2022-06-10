"""
WaterNetworkAnalysis
Set of tools for input preparation for conserved water search from MD trajectories (gromacs, amber) and their analysis
"""
import sys

from setuptools import find_packages, setup

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest-runner"] if needs_pytest else []

try:
    with open("README.rst") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name="WaterNetworkAnalysis",
    # version=0.0.1,
    author="Domagoj Fijan, Jelena Tosovic, Marko Jukic, Urban Bren",
    author_email="jecat_90@live.com",
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/x-rst",
    license="License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    keywords=("simulation analysis molecular dynamics biosimulation conserved water "),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    url="https://github.com/JecaTosovic/WaterNetworkAnalysis",
    download_url="https://pypi.org/project/WaterNetworkAnalysis/",
    project_urls={
        "Homepage": "https://github.com/JecaTosovic/WaterNetworkAnalysis",
        "Documentation": "https://WaterNetworkAnalysis.readthedocs.io/",
        "Source Code": "https://github.com/JecaTosovic/WaterNetworkAnalysis",
        "Issue Tracker": "https://github.com/JecaTosovic/WaterNetworkAnalysis/issues",
    },
    python_requires=">=3.6",
    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),
    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,
    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,
    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    install_requires=[
        "conservedwatersearch",
        "hdbscan>=0.8.28",
        "matplotlib",
        "MDAnalysis",
        "numpy",
        "scikit-learn",
        "nglview",
        "wget",
    ],  # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    zip_safe=False,
)
