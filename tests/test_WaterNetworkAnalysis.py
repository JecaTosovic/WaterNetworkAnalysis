"""
Unit and regression test for the WaterNetworkAnalysis package.
"""

import os

import MDAnalysis as mda
import numpy as np
import numpy.testing as npt

from WaterNetworkAnalysis import (
    align_trajectory,
    calculate_oxygen_density_map,
    extract_waters_from_trajectory,
    get_center_of_selection,
    get_selection_string_from_resnums,
    make_results_pdb_MDA,
)


def test_make_results_pdb_MDA():
    fname = "tests/data/MSR_HDBSCAN.dat"
    typefname = "tests/data/MSR_HDBSCAN_Type.dat"
    water_type = []
    waterO = []
    waterH1 = []
    waterH2 = []
    coords = np.loadtxt(fname)
    if coords.shape[1] == 3:
        for i in coords:
            waterO.append(i)
    else:
        Opos = coords[:, :3]
        H1 = coords[:, 3:6]
        H2 = coords[:, 6:9]
        for i, j, k in zip(Opos, H1, H2):
            waterO.append(i)
            waterH1.append(j)
            waterH2.append(k)
    types = np.loadtxt(typefname, dtype=str)
    for i in types:
        water_type.append(i)
    make_results_pdb_MDA(
        water_type=water_type,
        waterO=np.asarray(waterO),
        waterH1=np.asarray(waterH1),
        waterH2=np.asarray(waterH2),
        output_fname="test.pdb",
    )
    assert os.path.isfile("test.pdb")
    os.remove("test.pdb")


def test_get_center_of_selection():
    u = mda.fetch_mmtf("3t74")
    u.select_atoms("protein").write("test.pdb")
    sel = u.select_atoms("resid 123")
    cc = get_center_of_selection(get_selection_string_from_resnums([123]), "test.pdb")
    npt.assert_allclose(cc, sel.center(None))
    os.remove("test.pdb")


def test_calculate_oxygen_density_map():
    trjfname = "tests/data/testtrjgromacs.xtc"
    topfname = "tests/data/testtopgromacs.tpr"
    calculate_oxygen_density_map(
        get_center_of_selection("resname UBX", trjfname, topfname),
        trjfname,
        topfname,
    )
    os.remove("water.dx")
    os.remove("tests/data/.testtrjgromacs.xtc_offsets.npz")
    os.remove("tests/data/.testtrjgromacs.xtc_offsets.lock")


def test_extract_waters_from_trajectory():
    altrj = "tests/data/testalignedtrj.xtc"
    topof = "tests/data/testtopgromacs.tpr"
    extract_waters_from_trajectory(
        get_center_of_selection(
            get_selection_string_from_resnums(
                [
                    111,
                    112,
                    113,
                    122,
                    133,
                    138,
                    139,
                    142,
                    143,
                    157,
                    166,
                    167,
                    169,
                    170,
                    203,
                    231,
                    232,
                    238,
                ]
            ),
            altrj,
            topof,
        ),
        altrj,
        topof,
    )
    # os.remove("tests/data/testalignedtrj.xtc_offsets.npz")
    # os.remove("tests/data/testalignedtrj.xtc_offsets.lock")


def test_align_mda():
    align_trajectory(
        trajectory="tests/data/testtrjgromacs.xtc",
        output_trj_file="aligned_trajectory.xtc",
        align_target_file_name="align.pdb",
        align_mode="mda",
        align_target=0,
        align_selection="protein",
        topology="tests/data/testtopgromacs.tpr",
    )
    # os.remove("tests/data/testtrjgromacs_0.xtc")
    os.remove("aligned_trajectory.xtc")
    os.remove("align.pdb")


def test_align_mda_every():
    align_trajectory(
        trajectory="tests/data/testtrjgromacs.xtc",
        output_trj_file="aligned_trajectory.xtc",
        align_target_file_name="align.pdb",
        align_mode="mda",
        align_target=0,
        align_selection="protein",
        topology="tests/data/testtopgromacs.tpr",
        every=2,
    )
    # os.remove("tests/data/testtrjgromacs_0.xtc")
    os.remove("aligned_trajectory.xtc")
    os.remove("align.pdb")


def test_align_probis():
    align_trajectory(
        trajectory="tests/data/testtrjgromacs.xtc",
        output_trj_file="aligned_trajectory.xtc",
        align_target_file_name="align.pdb",
        align_mode="probis",
        align_target=0,
        align_selection="protein",
        topology="tests/data/testtopgromacs.tpr",
    )
    # os.remove("tests/data/testtrjgromacs_0.xtc")
    os.remove("aligned_trajectory.xtc")
    os.remove("probis")
    os.remove("align.pdb")


def test_density_map_units():
    trjfname = "tests/data/testalignedtrj.xtc"
    topfname = "tests/data/testtopgromacs.tpr"
    densgro = calculate_oxygen_density_map(
        get_center_of_selection("resname UBX", trjfname, topfname),
        trjfname,
        topfname,
    )
    trjpdbname = "tests/data/testalignedtrj.pdb"
    denspdb = calculate_oxygen_density_map(
        get_center_of_selection("resname UBX", trjpdbname),
        trjpdbname,
    )
    npt.assert_allclose(densgro.grid, denspdb.grid, atol=1e-3, rtol=1e-3)
