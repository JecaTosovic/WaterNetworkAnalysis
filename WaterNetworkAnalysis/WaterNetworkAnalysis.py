from __future__ import annotations

import os
import stat

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.density import DensityAnalysis, Density
from scipy.ndimage import gaussian_filter

from ConservedWaterSearch.utils import (
    get_orientations_from_positions,
    read_results,
)


def get_selection_string_from_resnums(
    resids: list[int], selection_type: str = "MDA"
) -> str:
    """Return selection string for given residue ids.

    Returns the selection command string for different programs based on
    amioacid residue IDs list given.

    Args:
        resids (list[int]): list of aminoacid residue ids.
        selection_type (str, optional): selection program language.
            Options:"MDA", "PYMOL", "PROBIS" and "NGLVIEW" Defaults to
            "MDA".

    Returns:
        str: selection command in form of a string

    Example::

        # list of resids
        resids = [8,12,143,144] # print PYMOL selection string
        get_selection_string_from_resnums(resids, selection_type =
        "PYMOL")
    """
    selection: str = ""
    for i in resids:
        if selection_type == "MDA":
            selection += "resnum " + str(i) + " or "
        if selection_type == "PYMOL":
            selection += str(i) + "+"
        if selection_type == "PROBIS":
            selection += str(i) + ","
        if selection_type == "NGLVIEW":
            selection += str(i) + " or "
    if selection_type == "MDA" or selection_type == "NGLVIEW":
        selection = selection[: len(selection) - 3]
    if selection_type == "PYMOL":
        selection = selection[: len(selection) - 1]
        selection = "resi " + selection
    if selection_type == "PROBIS":
        selection = selection[: len(selection) - 1]
        selection = (
            ' -motif1 "[:A and ('
            + selection
            + ')]" -motif2 "[:A and ('
            + selection
            + ')]" '
        )
    return selection


def get_center_of_selection(
    selection: str, trajectory: str, topology: str | None = None
) -> np.ndarray:
    """Compute centre of selection with MDAnalysis.

    Calculates coordinates in xyz of the centre of selection using
    MDAnalysis.

    Args:
        selection (str): selection string for MDAnalysis
        trajectory (str): trajectory filename - anything that is
            accepted by MDAnalysis should work
        topology (str | None, optional): topology filename. Defaults to
            None.

    Returns:
        np.ndarray:
            returns array that contains coordinates of center of selection

    Example::

        # find center of active site defined by residue ids
        resids = [8,12,143,144]
        get_center_of_selection(get_selection_string_from_resnums(resids))
    """
    if topology:
        u: mda.Universe = mda.Universe(topology, trajectory)
    else:
        u = mda.Universe(trajectory)
    s: mda.AtomGroup = u.select_atoms(selection)
    return s.center(None)


def calculate_oxygen_density_map(
    selection_center: np.ndarray,
    trajectory: str,
    topology: str | None = None,
    dist: float = 12.0,
    delta: float = 0.4,
    OW: str = "OW",
    SOL: str = "SOL",
    output_name: str = "water.dx",
) -> Density:
    """Generate oxygen density maps.

    Generate oxygen density maps using MDAnalysis.

    Args:
        selection_center (np.ndarray): center of selection
            around which waters will be selected.
        trajectory (str): trajectory filename.
        topology (str | None, optional): Topology filename if available.
            Defaults to None.
        dist (float, optional): distance around selection center inside
            which the oxygen will be selected. Defaults to 12.0.
        delta (float, optional): bin size for density map. Defaults to
            0.4 Angstroms.
        OW (str, optional): name of oxygens for selecting water oxygens.
            Defaults to "OW".
        SOL (str, optional): name of the solvent group for selecting
            solvent oxygens. Defaults to "SOL".
        output_name (str, optional): name of the output file, it should
            end with '.dx' . Defaults to "water.dx".

    Returns:
        Density:
            returns MDA Density object containing the density map

    Example::

        # Generate water oxygen density map near active site
        resids = [8,12,143,144] calculate_oxygen_density_map(
            get_center_of_selection(
                get_selection_string_from_resnums(resids)), trajectory =
                'trajectory.pdb'
            )
        )
    """
    if topology:
        u: mda.Universe = mda.Universe(topology, trajectory)
    else:
        u = mda.Universe(trajectory)
    ss: mda.AtomGroup = u.select_atoms(
        "name "
        + OW
        + " and resname "
        + SOL
        + " and point "
        + str(selection_center[0])
        + " "
        + str(selection_center[1])
        + " "
        + str(selection_center[2])
        + " "
        + str(dist),
        updating=True,
    )
    D: DensityAnalysis = DensityAnalysis(
        ss,
        delta=delta,
        gridcenter=selection_center,
        xdim=dist * 2,
        ydim=dist * 2,
        zdim=dist * 2,
    )
    D.run()
    # oxygen vdw radius
    vdw_radius = 1.52
    sigma = vdw_radius / delta
    probability_density = D.results.density.grid / D.results.density.grid.sum()
    smoothed_density = gaussian_filter(probability_density, sigma)
    smoothed_density = smoothed_density / np.max(smoothed_density)
    # Put smoothed density into a Density object
    smoothed_density_obj = Density(
        grid=smoothed_density,
        edges=D._edges,
        parameters={"isDensity": True},
    )
    smoothed_density_obj.export(output_name, type="double")
    return smoothed_density_obj


def make_results_pdb_MDA(
    water_type: list[str],
    waterO: np.ndarray,
    waterH1: np.ndarray,
    waterH2: np.ndarray,
    output_fname: str,
    protein_file: str | None = None,
    ligand_name: str | None = None,
    mode: str = "SOL",
) -> None:
    """Generate pdb file with clustering results.

    The water molecules determined by the clustering procedure are
    written in a pdb file which also contains protein and the ligand.
    Waters are labeled based on their hydrogen orientations (FCW for
    fully conserved, HCW for half conserved and WCW for weakly
    conserved). Uses MDAnalysis for construction of the pdb file.
    First 4 arguments of the function can be read from the results file
    by using cws.utils.read_results() or directly from the
    ``cws.water_clustering.WaterClustering`` class.

    Args:
        water_type (list[str]): List of water types.
        waterO (np.ndarray): numpy array containing coordinates
            of conserved waters' oxygens.
        waterH1 (np.ndarray): numpy array containing coordinates
            of conserved waters' first hydrogen
        waterH2 (np.ndarray): numpy array containing coordinates
            of conserved waters' second hydrogen
        output_fname (str): name of the output pdb file. Must end in '.pdb'.
        protein_file (str | None, optional): file which contains protein
            and ligand. It should be aligned in the same way as the
            trajectory used for calculation of conserved waters. If None no
            protein is saved. Defaults to None.
        ligand_name (str | None, optional): residue name of the ligand.
            If None no ligand is saved. Defaults to None.
        mode (str, optional): mode in which conserved waters will be
            saved. Options:

                - "SOL" - default mode. Saves water molecules as SOL so that
                  visualisation softwares can recognise them as waters. No
                  distinction is made between different types of conserved
                  waters.
                - "cathegorise" - cathegorises the waters according to
                  hydrogen orienation into fully conserved (FCW),
                  half-coserved (HCW) and weakly conserved (WCW). This mode
                  makes visualisers not able to recognise the waters as
                  water/sol but usefull for interpreting results.


    Example::

        # Generate pdb results file
        make_results_pdb_MDA(
            *cws.utils.read_results(),
            output_fname = 'results.pdb',
            mode = 'cathegorise'
        )
    """
    rescntr: int = 2
    # if protein file is not align_file then align the structure to the align file first??
    if protein_file is not None:
        us: mda.Universe = mda.Universe(protein_file)
        protein: mda.AtomGroup = us.select_atoms("protein")
        # here if ligand conformations have not been calculated perhaps try calculating them if possible?
        ligand: list = []
        if ligand_name:
            ligand.append(us.select_atoms("resname " + str(ligand_name)))
            # ne radi?
            # if (ligand[0]==None):
            #    print ("no ligand with given name "+str(ligand_name))
        rescntr += len(ligand) + len(protein.residues)
    # waters
    n_residues: int = len(water_type)
    n_atoms: int = n_residues * 3
    # create resindex list
    resindices: np.ndarray = np.repeat(range(n_residues), 3)
    # all water molecules belong to 1 segment
    segindices: list[int] = [0] * n_residues
    waters = mda.Universe.empty(
        n_atoms,
        n_residues=n_residues,
        atom_resindex=resindices,
        residue_segindex=segindices,
        trajectory=True,
    )  # necessary for adding coordintes
    waters.add_TopologyAttr("name", ["O", "H1", "H2"] * n_residues)
    waters.add_TopologyAttr("type", ["O", "H", "H"] * n_residues)
    # waters.add_TopologyAttr('segid', ['SOL'])
    coordinates: list[list[float]] = []
    rns: list[str] = []
    resid: list[int] = []
    looper = zip(waterO, waterH1, waterH2, water_type)
    for Opos, H1pos, H2pos, tip in looper:
        h2o = np.array(
            [
                [Opos[0], Opos[1], Opos[2]],  # oxygen
                [H1pos[0], H1pos[1], H1pos[2]],  # hydrogen
                [H2pos[0], H2pos[1], H2pos[2]],  # hydrogen
            ]
        )
        coordinates.extend(h2o)
        if mode == "cathegorise":
            if tip == "FCW":
                rns.append("FCW")
            if tip == "HCW":
                rns.append("HCW")
            if tip == "WCW":
                rns.append("WCW")
            if tip == "onlyO":
                rns.append("onlyO")
        elif mode == "SOL":
            rns.append("SOL")
        resid.append(rescntr)
        rescntr += 1
    waters.add_TopologyAttr("resname", rns)
    waters.add_TopologyAttr("resid", resid)
    coord_array = np.array(coordinates)
    waters.atoms.positions = coord_array
    bonds: list = []
    for o in range(0, n_atoms, 3):
        bonds.extend([(o, o + 1), (o, o + 2)])
    waters.add_TopologyAttr("bonds", bonds)
    if protein_file is not None:
        newU = mda.Merge(
            protein.atoms,
            *ligand,
            waters.atoms,
        )
    else:
        newU = mda.Merge(waters.atoms)
    newU.dimensions = 3
    if output_fname is None:
        if mode == "SOL":
            output_fname = "final_results_SOL.pdb"
        if mode == "cathegorise":
            output_fname = "final_results_cathegorise.pdb"
    newU.atoms.write(output_fname)


def extract_waters_from_trajectory(
    selection_center: np.ndarray,
    trajectory: str,
    topology: str | None = None,
    dist: float = 12.0,
    SOL: str = "SOL",
    OW: str = "OW",
    save_file: str | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract waters for clustering analysis.

    Calculates water (oxygen and hydrogen) coordinates for all the
    waters in the aligned trajectory using MDAnalysis for further use in water
    clustering. The trajectory should be aligned previously.

    Args:
        selection_center (np.ndarray): coordinates of selection
            center around which waters will be selected.
        trajectory (str): Trajectory file name.
        topology (str | None, optional): Topology file name. Defaults to None.
        dist (float, optional): Distance around the center of selection
            inside which water molecules will be sampled. Defaults to 12.0.
        SOL (str, optional): Residue name for waters. Defaults to "SOL".
        OW (str, optional): Name of the oxygen atom in water molecules.
            Defaults to "OW".
        save_file (str | None, optional): File to which coordinates will
            be saved. If none doesn't save to a file. Defaults to None.

    Returns:
        tuple[np.ndarray, np.ndarray]:
            returns xyz numpy arrays that contain coordinates of oxygens,
            and combined array of hydrogen 1 and hydrogen 2

    Example::

        # Generate water coordinates for clustering analysis
        resids = [8,12,143,144]
        coordO, coordH = extract_waters_from_trajectory(
            get_center_of_selection(get_selection_string_from_resnums(resids)),
            trajectory = 'trajectory.xtc',
            topology = 'topology.tpr'
        )
    """
    if topology:
        u = mda.Universe(topology, trajectory)
    else:
        u = mda.Universe(trajectory)
    coordsH = []
    coordsO = []
    waters = u.select_atoms(f"resname {SOL}")
    if len(waters.bonds) == 0:
        waters.guess_bonds()
    # loop over
    for nn, k in enumerate(u.trajectory):
        Os = u.select_atoms(
            "name "
            + str(OW)
            + " and resname "
            + str(SOL)
            + " and point "
            + str(selection_center[0])
            + " "
            + str(selection_center[1])
            + " "
            + str(selection_center[2])
            + " "
            + str(dist)
        )
        for i, j in zip(Os.positions, Os):
            Hs = j.bonded_atoms
            if len(Hs) != 2:
                raise Exception(
                    f"Water {j} in snapshot {i} has too many hydrogens ({len(Hs)})."
                )
            for hydrogen_positions in Hs.positions:
                coordsH.append(hydrogen_positions)
            coordsO.append(i)
    Odata: np.ndarray = np.asarray(coordsO)
    coordsH = np.asarray(coordsH)
    Opos, H1, H2 = get_orientations_from_positions(Odata, coordsH)
    # SAVEs full XYZ coordinates, not O coordinates and h orientations!!!!!
    if save_file is not None:
        np.savetxt(save_file, np.c_[Opos, H1, H2])
    return Odata, coordsH


def __align_mda(
    unaligned_trj_file: str,
    pdb_to_align_to: str,
    output_trj_file: str,
    topology: str | None = None,
    align_selection: str = "protein",
) -> None:
    """Alignment via MDAnalysis.

    This function is not ment to be used directly. Use align_trajectory instead.

    Args:
        unaligned_trj_file (str): trajectory file containing unaligned
            trajectory.
        pdb_to_align_to (str): pdb file that containes the reference
            state to align to.
        output_trj_file (str): output file trajectory will be saved to.
        topology (str | None, optional): File containing topology.
            Defaults to None.
        align_selection (str, optional): selection to align to. Defaults
            to "protein".
    """
    if topology:
        mob: mda.Universe = mda.Universe(topology, unaligned_trj_file)
    else:
        mob = mda.Universe(unaligned_trj_file)
    ref: mda.Universe = mda.Universe(pdb_to_align_to)
    mda.analysis.align.AlignTraj(
        mob, ref, select=align_selection, filename=output_trj_file
    ).run()


def __align_probis(
    unaligned_trj_file: str,
    pdb_to_align_to: str,
    output_trj_file: str,
    topology: str | None = None,
    align_selection: str = "protein",
    probis_exec: str | None = None,
) -> None:
    """Alignment function that uses probis.

    This function is not ment to be used directly. Use align_trajectory
    instead.Uses probis algorithm for alignment. See
    `this link <http://insilab.org/probis-algorithm/>`__ for more info.

    Args:
        unaligned_trj_file (str): trajectory file containing unaligned
            trajectory.
        pdb_to_align_to (str): pdb file that containes the reference
            state to align to.
        output_trj_file (str): output file trajectory will be saved to.
        topology (str | None, optional): File containing topology.
            Defaults to None.
        align_selection (str, optional): selection to align to. Defaults
            to "protein".
        probis_exec (str | None, optional): Executable to run probis
            algorithm. If None it is downloaded from the internet.
            Defaults to None.
    """
    import wget

    if probis_exec is None:
        pwd: str = os.path.abspath(os.getcwd())
        probis_web: str = "http://insilab.org/files/probis-algorithm/probis"
        if not (os.path.isfile(pwd + "/probis")):
            wget.download(probis_web, pwd + "/probis")
        probis_exec = "probis"
        st: os.stat_result = os.stat(probis_exec)
        os.chmod(probis_exec, st.st_mode | stat.S_IEXEC)
    if align_selection == "protein":
        sele: str = ""
    else:
        sele = align_selection
    if topology:
        u: mda.Universe = mda.Universe(topology, unaligned_trj_file)
    else:
        u = mda.Universe(unaligned_trj_file)
    mobile_pdb: str = "mobile.pdb"
    u.select_atoms(align_selection).segments.segids = "A"
    u.add_TopologyAttr("chainIDs")
    u.select_atoms(align_selection).chainIDs = "A"
    forwrite = u.atoms
    with mda.Writer(output_trj_file, multiframe=True) as W:
        for snap in u.trajectory:
            with mda.Writer(mobile_pdb) as mw:
                mw.write(forwrite)
            # string to run
            runstring: str = (
                probis_exec
                + " -compare -super "
                + sele
                + " -f1 "
                + str(pdb_to_align_to)
                + " -c1 A -f2 "
                + str(mobile_pdb)
                + " -c2 A"
            )
            if probis_exec == "probis":
                runstring = "./" + runstring
            # run probis
            print(runstring)
            os.system(runstring)
            # read in rota.pdb files and remove all other files
            newal: mda.Universe = mda.Universe(
                pdb_to_align_to[:4] + "A_" + mobile_pdb[:4] + "A.0.rota.pdb"
            )
            # remove aligned pdb files
            removestring: str = "rm *.rota.pdb probis*.tmp " + mobile_pdb
            os.system(removestring)
            # add cell dimensions and other missing info?
            newal.dimensions = snap.dimensions
            # append snapshot to new trajectory
            selnew = newal.atoms
            W.write(selnew)


def align_trajectory(
    trajectory: str,
    output_trj_file: str,
    align_target_file_name: str,
    topology: str | None = None,
    every: int = 1,
    align_mode: str = "mda",
    align_target: int | None = -1,
    align_selection: str = "protein",
    probis_exec: str | None = None,
) -> None:
    """Align the trajectory.

    Before running water clustering for identification of conserved
    water molecules the trajectory should be aligned first. Alignment
    can be done via MDAnalysis or using the probis algorithm. Whole
    protein is aligned by default. To select the align reference state
    either select an integer for ``align_target`` and specify a file name to
    which the align target will be saved to with ``align_target_file_name``
    OR set align_target to ``None`` and ``align_target_file_name`` will be read
    and used as align target.

    The trajectory or topology should contain information on bond topology
    for alignment. Supported topology file types:

        DATA
        DMS
        GSD
        MMTF
        MOL2
        PARMED
        PDB
        ENT
        PSF
        TOP
        PRMTOP
        PARM7
        TPR
        TXYZ
        ARC
        XML
        XPDB

    Alternatively the whole trajectory can be provided in some of the above
    given file types as well.

    Args:
        trajectory (str): File name containing unaligned trajectory.
        output_trj_file (str): output file name for aligned trajectory.
        align_target_file_name (str): File name for saving the align
            target (usually pdb) if ``align_target`` is int. If align target is
            ``None``, the align target will be read from this file instead.
        topology (str | None, optional): Topology file name. Defaults to
            ``None``.
        every (int, optional): Take every ``every`` snapshot instead
            of taking all the snapshots (every = 1) for alignment. Defaults
            to 1.
        align_mode (str, optional): Align algorithm to use. "mda"
            uses MDAnalysis while "probis" uses the probis algorithm.
            Defaults to "mda".
        align_target (int | None, optional): Align target. If ``None`` the
            align target is read from the ``align_target_file_name``. If a
            number is given uses the given snapshot of the trajectory as
            the align target. If -1 uses the last snapshot. Defaults to
            -1.
        align_selection (str, optional): Selection to align to. Defaults
            to "protein".
        probis_exec (str | None, optional): location of probis
            executable if probis is used. If ``None`` it is downloaded from the
            internet. Defaults to ``None``.

    Example::

        # align the trajectory and save to a file
        align_trajectory(
            trajectory="trajectory.xtc",
            output_trj_file="aligned_trajectory.xtc",
            align_target_file_name='aligned.pdb', align_mode="mda",
            align_target=0, align_selection="protein",
            topology="topology.tpr",
        )
    """
    import MDAnalysis.transformations as trans

    if topology is not None:
        mob = mda.Universe(topology, trajectory)
        ref: mda.Universe = mda.Universe(topology, trajectory)
    else:
        mob: mda.Universe = mda.Universe(trajectory)
        ref = mda.Universe(trajectory)
    # center the box to center of mass
    if topology is None:
        if not (
            trajectory.upper().endswith(
                (
                    "DATA",
                    "DMS",
                    "GSD",
                    "MMTF",
                    "MOL2",
                    "PARMED",
                    "PDB",
                    "ENT",
                    "PSF",
                    "TOP",
                    "PRMTOP",
                    "PARM7",
                    "TPR",
                    "TXYZ",
                    "ARC",
                    "XML",
                    "XPDB",
                    "PDB",
                )
            )
        ):
            raise Exception(
                "unsupported trajectory file format. bond topology information is needed for alignment. Consider providing a topology file."
            )
    elif not (
        topology.upper().endswith(
            (
                "DATA",
                "DMS",
                "GSD",
                "MMTF",
                "MOL2",
                "PARMED",
                "PDB",
                "ENT",
                "PSF",
                "TOP",
                "PRMTOP",
                "PARM7",
                "TPR",
                "TXYZ",
                "ARC",
                "XML",
                "XPDB",
                "PDB",
            )
        )
    ):
        raise Exception(
            "unsupported topology file type. Bond information is needed for alignment."
        )
    ref.select_atoms(align_selection).segments.segids = "A"
    ref.add_TopologyAttr("chainIDs")
    ref.select_atoms(align_selection).chainIDs = "A"
    if type(align_target) == int:
        wr = ref.atoms
        wr.write(align_target_file_name, frames=ref.trajectory[[align_target]])
    if every > 1:
        unaligned_trj_file: str = (
            trajectory[: trajectory.rfind(".") - 1]
            + f"every_{every}"
            + trajectory[trajectory.rfind(".") :]
        )
        with mda.Writer(unaligned_trj_file, multiframe=True) as W:
            for i in mob.trajectory[::every]:
                W.write(mob.atoms)
    elif every == 1:
        unaligned_trj_file: str = trajectory
    else:
        raise Exception("every must be a positive integer")
    if align_mode == "probis":
        __align_probis(
            unaligned_trj_file=unaligned_trj_file,
            pdb_to_align_to=align_target_file_name,
            output_trj_file=output_trj_file,
            topology=topology,
            align_selection=align_selection,
            probis_exec=probis_exec,
        )
    elif align_mode == "mda":
        __align_mda(
            unaligned_trj_file=unaligned_trj_file,
            pdb_to_align_to=align_target_file_name,
            output_trj_file=output_trj_file,
            topology=topology,
            align_selection=align_selection,
        )
    else:
        raise Exception("wrong alignemnt mode, must be probis or mda")


def align_and_extract_waters(
    center_for_water_selection: np.ndarray,
    trajectory: str,
    aligned_trajectory_filename: str,
    align_target_file_name: str,
    topology: str | None = None,
    every: int = 1,
    align_mode: str = "mda",
    align_target: int | None = -1,
    align_selection: str = "protein",
    probis_exec: str | None = None,
    dist: float = 12.0,
    SOL: str = "SOL",
    OW: str = "OW",
    HW1: str = "HW1",
    HW2: str = "HW2",
) -> tuple[np.ndarray]:
    """Align and extracts waters from trajectory.

    Aligns the trajectory first and then extracts water molecules for
    further water clustering analysis. If trajectory has already been
    aligned, one can use :py:meth:`extract_waters_from_trajectory` to
    extract the water molecules for water clustering analysis.

    Args:
        center_for_water_selection (np.ndarray): Coordiantes
            around which all water molecules inside a radius ``dist``
            will be seleceted for water clustering analysis.
        trajectory (str): File name of the trajectory from which waters
            will be extracted.
        aligned_trajectory_filename (str): File name to which aligned
            trajectory will be saved.
        align_target_file_name (str): File name for saving the align
            target (usually pdb) if ``align_target`` is ``int``. If
            align target is ``None``, the align target will be read
            from this file instead!
        topology (str | None, optional): Topology file name. Defaults
            to ``None``.
        every (int, optional): Take every ``every`` snapshot instead of
            taking all the snapshots (every = 1) for alignment.
            Defaults to 1.
        align_mode (str, optional): Align algorithm to use. "mda" uses
            MDAnalysis while "probis" uses the probis algorithm.
            Defaults to "mda".
        align_target (int | None, optional): Align target. If ``None``
            the align target is read from the align_target_file_name.
            If a number is given uses the given snapshot of the
            trajectory as the align target. If -1 uses the last
            snapshot. Defaults to -1.
        align_selection (str, optional): Selection to align to. Defaults
            to "protein".
        probis_exec (str | None, optional): location of probis
            executable if probis is used. If None it is downloaded
            from the internet. Defaults to None.
        dist (float, optional): Radius around
            ``center_for_water_selection`` to be used for extraction of
            water molecules. Defaults to 12.0.
        SOL (str, optional): Residue name for waters. Defaults to "SOL".
        OW (str, optional): Name of the oxygen atom in water molecules. Defaults
             to "OW".
        HW1 (str, optional): Name of hydrogen 1 atom in water
            molecules. Defaults to "HW1".
        HW2 (str, optional): Name of hydrogen 2 atom in water molecules.
            Defaults to "HW2".
        save_file (str | None, optional): File to which coordinates
            will be saved. If none doesn't save to a file. Defaults to None.

    Returns:
        tuple[np.ndarray, np.ndarray]:
            Returns coordinates of oxygen atoms, first
            hydrogen atom and second hydrogen atom in three seperate numpy
            arrays. Each row in each array makes up coordinates of a single
            water molecule.

    Example::

        # Generate water coordinates for clustering analysis from unaligned trajectory
        resids = [8,12,143,144] align_and_extract_waters(
            get_center_of_selection(get_selection_string_from_resnums(resids)),
            trajectory = 'trajectory.xtc', aligned_trajectory_filename =
            'aligned_trj.xtc', align_target_file_name = 'aligned.pdb',
            topology = 'topology.tpr', every = 1, align_mode = "mda",
            align_target= 0, align_selection = "protein", dist = 10.0,
        )
    """
    align_trajectory(
        trajectory=trajectory,
        output_trj_file=aligned_trajectory_filename,
        align_target_file_name=align_target_file_name,
        topology=topology,
        every=every,
        align_mode=align_mode,
        align_target=align_target,
        align_selection=align_selection,
        probis_exec=probis_exec,
    )
    return extract_waters_from_trajectory(
        center_for_water_selection,
        trajectory=aligned_trajectory_filename,
        dist=dist,
        topology=topology,
        SOL=SOL,
        OW=OW,
        HW1=HW1,
        HW2=HW2,
    )


def read_results_and_make_pdb(
    fname: str = "Clustering_results.dat",
    typefname: str = "Type_Clustering_results.dat",
    protein_file: str = "aligned.pdb",
    ligand_name: str | None = None,
    output_fname: str | None = None,
    mode: str = "SOL",
) -> None:
    """Read results from files and generate a pdb file.

    Args:
        fname (str, optional): File name with clustering results -
            coordinates of water molecules. Defaults to "Clustering_results.dat".
        typefname (str, optional): File name which contains the types of
            each water molecule. Defaults to "Type_Clustering_results.dat".
        protein_file (str, optional): File name which contains the
            reference structure trajectory has been aligned to. Defaults
            to "aligned.pdb".
        ligand_name (str | None, optional): Residue name for the ligand.
            If none is given, no ligand is extracted and
            visualised/saved in the pdb file. Defaults to None.
        output_fname (str | None, optional): Name of the output file
            (pdb prefered). Defaults to None.
        mode (str, optional): mode in which conserved waters will
            besaved. Options:

                - "SOL" - default mode. Saves water molecules as SOL so that
                  visualisation softwares can recognise them as waters. No
                  distinction is made between different types of conserved
                  waters.
                - "cathegorise" - cathegorises the waters according to
                  hydrogen orienation into fully conserved (FCW),
                  half-coserved (HCW) and weakly conserved (WCW). This mode
                  makes visualisers not able to recognise the waters as
                  water/sol but usefull for interpreting results.

    Example::

        # Generate pdb results file
        read_results_and_make_pdb(
            fname = 'Results.dat',
            typefname = 'TypeResults.dat',
            protein_file = 'aligned_protein.pdb',
            ligand_name = 'UBY',
            output_fname = 'results.pdb',
            mode = 'cathegorise',
        )
    """
    water_type, waterO, waterH1, waterH2 = read_results(fname, typefname)
    make_results_pdb_MDA(
        water_type=water_type,
        waterO=waterO,
        waterH1=waterH1,
        waterH2=waterH2,
        protein_file=protein_file,
        ligand_name=ligand_name,
        output_fname=output_fname,
        mode=mode,
    )
