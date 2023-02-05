"""
WaterNetworkAnalysis
Set of tools for input preparation for conserved water search from MD trajectories (gromacs, amber) and their analysis
"""
import sys

from setuptools import find_packages, setup

short_description = __doc__.split("\n")

needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest-runner"] if needs_pytest else []

try:
    with open("README.rst") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])


setup(
    name="WaterNetworkAnalysis",
    version="0.1.0",
    author="Domagoj Fijan, Jelena Tosovic, Marko Jukic, Urban Bren",
    author_email="jecat_90@live.com",
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/x-rst",
    license="License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    keywords=("simulation analysis molecular dynamics biosimulation conserved water "),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    url="https://github.com/JecaTosovic/WaterNetworkAnalysis",
    download_url="https://pypi.org/project/WaterNetworkAnalysis/",
    project_urls={
        "Homepage": "https://github.com/JecaTosovic/WaterNetworkAnalysis",
        "Documentation": "https://WaterNetworkAnalysis.readthedocs.io/",
        "Source Code": "https://github.com/JecaTosovic/WaterNetworkAnalysis",
        "Issue Tracker": "https://github.com/JecaTosovic/WaterNetworkAnalysis/issues",
    },
    python_requires=">=3.8",
    packages=find_packages(),
    include_package_data=True,
    setup_requires=[] + pytest_runner,
    install_requires=[
        "ConservedWaterSearch",
        "MDAnalysis",
        "numpy<1.24",
        "wget",
    ],
    platforms=["Linux", "Mac OS-X", "Unix", "Windows"],
    zip_safe=False,
)
