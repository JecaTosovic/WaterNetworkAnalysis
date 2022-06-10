"""
WaterNetowrkAnalysis.py
Module for analysis of surface water networks near protein active sites.

Handles the primary functions
"""

import os
import shutil
import stat
import sys
from collections.abc import Iterable
from typing import List

import hdbscan
import matplotlib.pyplot as plt
import MDAnalysis as mda
import MDAnalysis.analysis
import MDAnalysis.analysis.align
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import OPTICS, KMeans, cluster_optics_dbscan, cluster_optics_xi

import conservedwatersearch
from conservedwatersearch.hydrogen_orientation import *
from conservedwatersearch.utils import *


class WaterNetworkAnalysis:
    """
    Main class for surface water analysis. Contains functions for surface water analysis workflow which usually consists of a) finding waters by first alligning a given trajectory and selecting relevent waters; b) preforming anaylisis such as clustering of water networks of selected waters.

    Parameters
    ----------
    aminoacids_in_activesite : list[int], Optional
        list of ints with aminoacid residue ids which belong to active site residues.

    Attributes
    -------
    """

    randomcntr = 0
    plotend = 0
    plotter = 0
    debugH = 0
    debugO = 0
    verbose = 0
    other_waters_minsamp_pct = 0.15
    # controls percentage of minsamp devation from number of elements in cluster of oxygens for other waters
    Hydxi = 0.05
    # xi for optics for clustering of hydrogen orientations
    kmeansinertiacutoff = 0.4
    kmeansangcutoff = 120
    conserved_optics_angdiff_cutoff = 5
    conserved_optics_angstd_cutoff = 17
    half_conserved_angle_limit = 15
    wiggly_conserved_std_limit = 20
    wiggly_explained = 0.7
    numbpct = 0.80
    active_site_center = 0
    lower_minsamp = 0.75
    njobs = 1
    pdbfile = "1.pdb"
    maxsamp = None
    minsamp = 50
    analyse_hydrogen = True
    delete_after_find = True
    xi = 0.01
    X = []
    clust = 0
    space = 0
    reachability = 0
    labels = 0
    eta = 0
    clustcoordO = []
    clustnames = []
    conservedO = []
    conservedH1 = []
    conservedH2 = []
    all_waters = []
    waterO = []
    waterH1 = []
    waterH2 = []
    water_type = []
    kmeans = KMeans(n_clusters=2)
    waters = []
    aminokis_u_am = []
    activesite_aminoacids = []
    activesite_aminoacid_ids = []
    H = []
    H1 = []
    H2 = []
    v1 = []
    v2 = []
    origin_water = 0
    isDBSCAN = False
    isHDBSCAN = False
    ligand_conf = None
    align_file = None
    maxsamp = None
    clustering_type = None
    clustering_options = None
    clustering_algorithm = None
    nsnaps = None
    target_oxygen_cluster_size = None
    numbpct_hyd_orient_analysis = 0.85

    def __init__(self, aminoacids_in_activesite: List[int] = None):
        if aminoacids_in_activesite:
            self.set_aminoacids_in_activesite(aminoacids_in_activesite)

    def set_aminoacids_in_activesite(self, aminoacids_in_activesite: List[int]):
        """
        Sets aminoacids in active site for future water extraction and other analysis.

        Parameters
        ----------
        aminoacids_in_activesite : List[int]
            list of ints with aminoacid residue ids which belong to active site residues.
        """
        self.activesite_aminoacid_ids = aminoacids_in_activesite
        self.activesite_aminoacids_mda = ""
        self.activesite_aminoacids = ""
        for i in aminoacids_in_activesite:
            self.activesite_aminoacids_mda += "resnum " + str(i) + " or "
            self.activesite_aminoacids += str(i) + "+"
        self.activesite_aminoacids_mda = self.activesite_aminoacids_mda[
            : len(self.activesite_aminoacids_mda) - 3
        ]
        self.activesite_aminoacids = self.activesite_aminoacids[
            : len(self.activesite_aminoacids) - 1
        ]

    def setdebuggers(self):
        """
        Sets the debugger options to avoid clashes.
        """
        np.set_printoptions(precision=4)
        if self.plotend == 1:
            self.debugO = 1
            self.debugH = 1
        if (self.debugO == 1 or self.debugH == 1) and self.verbose < 1:
            self.verbose = 1
        if self.debugO == 2 or self.debugH == 2:
            self.verbose = 2

    def calculate_oxygen_density_map(
        self,
        trajectory: str,
        topology: str = None,
        dist: float = 12.0,
        OW: str = "OW",
        SOL: str = "SOL",
        output_name: str = "water.dx",
    ):
        """
        Generates (water)oxygen density map.

        Parameters
        ----------
        trajectory : str
            trajectory filename - make sure its aligned!
        topology : str, optional
            topology filename - tpr for gromacs prefered, by default None
        dist : float, optional
            distance from active site center to select waters in, by default 12.0
        OW : str, optional
            name of oxygens for selecting water oxygens, by default "OW"
        SOL : str, optional
            name solvent group for selecting solvent oxygens, by default "SOL"
        output_name : str, optional
            name of outputfile - it should be .dx, by default "water.dx"
        """
        from MDAnalysis.analysis.density import DensityAnalysis

        if topology:
            u = mda.Universe(topology, trajectory)
        else:
            u = mda.Universe(trajectory)
        s = u.select_atoms(self.activesite_aminoacids_mda)
        cc = s.center(None)
        ss = u.select_atoms(
            "name "
            + OW
            + " and resname "
            + SOL
            + " and point "
            + str(cc[0])
            + " "
            + str(cc[1])
            + " "
            + str(cc[2])
            + " "
            + str(dist)
        )
        D = DensityAnalysis(ss)
        D.run()
        D.density.convert_density("A^{-3}")
        D.density.export(output_name, type="double")

    def visualise_pymol(
        self,
        pdbfile: str = "final_results_SOL.pdb",
        pse_file: str = "clustering_session.pse",
        create_pse_only: bool = True,
        crystal_waters: str = None,
        density_map: str = None,
        ligand_resname: str = None,
    ):
        """
        Starts a pymol visualisation instance by reading clustering result made by make_results_pdb and saves it to pse_file

        Parameters
        ----------
        pdbfile : str, default : "final_results_SOL.pdb"
            input pdb file which will be used for visualisation; it should contain the protein, ligand/lignads and surface waters. The default is pdb file generated from make_results_pdb.
        pse_file : str, default: "clustering_session.pse"
            file in which to save the clustering session generated by this function.
        create_pse_only, bool, default: True
            only creates pse file without starting visualisation - useful if your pymol installation is not stable in linux but you still want pse file created.
        crystal_waters, str, default: None
            PDB database ID from which the waters will be added to the pymol session, if you define ligand_resname crystal waters will be selected around 10 A from the resname, default None - don't add crystal waters
        density_map, str, default: None
            .dx file from which the volume representation of density will be plotted. Make sure that you produce this from same aligned trajectory as used for analysis, default None - don't add density volumes
        ligand_resname, str, default: None
            if you define ligand_resname crystal waters will be selected around 10 A from the resname, default None - all waters will be selected

        """
        import pymol
        from pymol import cmd

        if not (create_pse_only):
            pymol.finish_launching(["pymol", "-q"])
        cmd.load(pdbfile)
        cmd.hide("everything")
        # polymer for surface def
        tmpObj = cmd.get_unused_name("_tmp")
        cmd.create(tmpObj, "( all ) and polymer", zoom=0)
        # aminoacids in active site
        aminokis_u_am = cmd.get_unused_name("ak_u_am_")
        cmd.select(
            aminokis_u_am, " resi " + self.activesite_aminoacids + " and polymer "
        )
        cmd.show("licorice", aminokis_u_am)
        cmd.color("magenta", aminokis_u_am)
        # pseudoatom in active site center
        active_site_center = cmd.get_unused_name("centaraktivnogmesta_")
        cmd.pseudoatom(active_site_center, aminokis_u_am)
        cmd.hide(representation="everything", selection=active_site_center)
        # cmd.show("sphere",active_site_center)
        # cmd.set ("sphere_scale",0.1,active_site_center)
        # protein surface
        protein = cmd.get_unused_name("samo_protein_")
        cmd.select(protein, "polymer")
        povrsina = cmd.get_unused_name("povrsina_protein_")
        cmd.create(povrsina, protein)
        cmd.show("surface", povrsina)
        cmd.color("gray70", povrsina)
        cmd.set("transparency", 0.5, povrsina)
        # ligand representation
        ligand = cmd.get_unused_name("ligand_")
        cmd.select(ligand, "organic")
        cmd.show("licorice", ligand)
        # add water representations
        if "cathegorise" in pdbfile:
            # conserved waters
            conserved = cmd.get_unused_name("conserved_")
            cmd.select(conserved, "resname FCW")
            cmd.show("lines", conserved)
            cmd.color("firebrick", conserved)
            # half conserved waters
            half_conserved = cmd.get_unused_name("half_con_")
            cmd.select(half_conserved, "resname HCW")
            cmd.show("lines", half_conserved)
            cmd.color("skyblue", half_conserved)
            # semi conserved waters
            semi_conserved = cmd.get_unused_name("wig_con_")
            cmd.select(semi_conserved, "resname WCW")
            cmd.show("lines", semi_conserved)
            cmd.color("limegreen", semi_conserved)
            # all waters
            waters = cmd.get_unused_name("waters_")
            cmd.select(
                waters, semi_conserved + " or " + half_conserved + " or " + conserved
            )
        else:
            waters = cmd.get_unused_name("waters_")
            cmd.select(waters, "SOL")
            cmd.show("lines", waters)
            # distances for hydrogen bonds
            cmd.distance(
                "polar_contacts",
                aminokis_u_am + " or " + waters + " or organic",
                "sol",
                mode=2,
            )
            cmd.hide("labels")
        # Add crystal waters
        if crystal_waters:
            cmd.fetch(crystal_waters)
            cmd.hide("everything", crystal_waters)
            cmd.align(
                "polymer and " + crystal_waters, "polymer and " + pdbfile.split(".")[0]
            )
            if ligand_resname:
                cmd.select(
                    "crystal_waters",
                    "("
                    + crystal_waters
                    + " and SOL) within 10 of resname "
                    + ligand_resname,
                )
            else:
                cmd.select("crystal_waters", crystal_waters + " and SOL")
            cmd.show("spheres", "crystal_waters")
            cmd.set("sphere_scale", "0.4", "crystal_waters")
        # add volume density visualisation
        if density_map:
            cmd.load(density_map)
            cmd.volume("water_density", density_map.split(".")[0])
        # reset camera
        cmd.reset()
        cmd.center(active_site_center)
        # save
        cmd.save(pse_file)
        cmd.reinitialize()

    def visualise_ngl(self, pdbfile: str = "final_results_SOL.pdb"):
        """
        Starts a nglview visualisation instance by reading clustering result made by make_results_pdb. Requires jupyter notebook.

        Parameters
        ----------
        pdbfile : str, default : "final_results_SOL.pdb"
            input pdb file which will be used for visualisation; it should contain the protein, ligand/lignads and surface waters. The default is pdb file generated from make_results_pdb.

        Returns
        -------
        view : nglview
            returns nglview instances which should be saved and used in single line to start visualisation.
        """
        import nglview as ngl

        view = ngl.show_mdanalysis(mda.Universe(pdbfile), gui=True)
        view.clear_representations()
        view.add_surface(selection="protein", opacity=0.5)
        view.add_representation("ball+stick", selection="not protein and not HC water")
        # view.add_representation('ball+stick', selection='not protein and not HC and not SC and not CON')
        # view.add_representation('ball+stick', selection='CON',color='red')
        # view.add_representation('ball+stick', selection='HC',color='blue')
        # view.add_representation('ball+stick', selection='SC',color='green')
        view.add_representation("ball+stick", selection="water", color="red")
        view.add_representation(
            "ball+stick",
            selection=self.activesite_aminoacids.replace("+", " or "),
            color="resname",
        )
        # can I add distance for hydrogen bonds here maybe as well?
        # view.add_representation("distance",)
        # view.add_representation("contact",)
        return view

    def make_results_pdb(
        self,
        protein_file: str = None,
        ligand_name: str = None,
        output_fname: str = None,
        mode: str = "SOL",
    ):
        """
        Makes a pdb file from results of surface water network clustering that contains protein, ligand and surface waters. The waters can be saved as "SOL" in pdb file or as the conserved ("CON"), half-conserved ("HC") and semi-conserved ("SC") depending on mode. The generated file can be used by visualise_ngl and visualise_pymol functions.

        Parameters
        ----------
        protein_file : str, default : self.align_file
            input pdb file which will be used for extraction of protein and ligand. It should be aligned to the aligned trajectory. Default is the conformation to which the trajectory is aligned.
        ligand_name : str, Optional, default : self.ligand_conf
            residue name of the ligand. used for extracting the ligand geometry from the pdb. If None ?? use the conformations extracted from ??
        output_fname : str, default : "final_results_$.pdb"
            Output pdb filename. Default name is based on mode - $ replaces the mode.
        mode : str, default : "SOL"
            Mode for water residue representation in pdb file. Possible choices:
                "SOL" - default mode. Saves water molecules as SOL so that visualisation softwares can recognise them as waters. No distinction is made between different types of conserved waters.
                "cathegorise" - cathegorises the waters according to hydrogen orienation into conserved (CON), half-coserved (HC) and semi-conserved (SC). This mode makes visualisers not able to recognise the waters as water/sol. use with care.

        See Also
        --------
        visualise_ngl, visualise_pymol
        """
        if protein_file is None:
            protein_file = self.align_file
        # if protein file is not self.align_file then align the structure to the align file first??
        us = mda.Universe(protein_file)
        protein = us.select_atoms("protein")
        # here if ligand conformations have not been calculated perhaps try calculating them if possible?
        ligand = []
        # ne radi?
        # print(self.ligand_conf)
        # if self.ligand_conf:
        if False:
            pass
        else:
            if ligand_name:
                ligand.append(us.select_atoms("resname " + str(ligand_name)))
                # ne radi?
                # if (ligand[0]==None):
                #    print ("no ligand with given name "+str(ligand_name))
        # waters
        n_residues = len(self.water_type)
        n_atoms = n_residues * 3
        # create resindex list
        resindices = np.repeat(range(n_residues), 3)
        # all water molecules belong to 1 segment
        segindices = [0] * n_residues
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
        coordinates = []
        rns = []
        resid = []
        rescntr = len(protein.residues) + len(ligand) + 2
        looper = zip(self.waterO, self.waterH1, self.waterH2, self.water_type)
        for O, H1, H2, tip in looper:
            h2o = np.array(
                [
                    [O[0], O[1], O[2]],  # oxygen
                    [H1[0], H1[1], H1[2]],  # hydrogen
                    [H2[0], H2[1], H2[2]],
                ]
            )  # hydrogen
            coordinates.extend(h2o)
            if mode == "cathegorise":
                if tip == "FCW":
                    rns.append("FCW")
                if tip == "HCW":
                    rns.append("HCW")
                if tip == "WCW":
                    rns.append("WCW")
                if tip == "O_clust":
                    rns.append("OXY")
            elif mode == "SOL":
                rns.append("SOL")
            resid.append(rescntr)
            rescntr += 1
        waters.add_TopologyAttr("resname", rns)
        waters.add_TopologyAttr("resid", resid)
        coord_array = np.array(coordinates)
        assert coord_array.shape == (n_atoms, 3)
        waters.atoms.positions = coord_array
        bonds = []
        for o in range(0, n_atoms, 3):
            bonds.extend([(o, o + 1), (o, o + 2)])
        waters.add_TopologyAttr("bonds", bonds)
        newU = mda.Merge(
            protein.atoms,
            *ligand,
            waters.atoms,
        )
        newU.dimensions = us.dimensions
        if output_fname is None:
            if mode == "SOL":
                output_fname = "final_results_SOL.pdb"
            if mode == "cathegorise":
                output_fname = "final_results_cathegorise.pdb"
        newU.atoms.write(output_fname)

    def add_visualisation_waters(self, waters):
        """Adds waters to variables which will be used for visuallising waters"""
        for i in waters:
            self.waterO.append(i[0])
            if len(i) > 2:
                self.waterH1.append(i[1])
                self.waterH2.append(i[2])
            self.water_type.append(i[-1])

    def extract_waters_from_trajectory(
        self,
        traj: str = "aligned_trajectory.xtc",
        topology: str = None,
        dist: float = 12.0,
        SOL: str = "SOL",
        OW: str = "OW",
        HW1: str = "HW1",
        HW2: str = "HW2",
    ):
        """
        Function that extracts waters from a trajectory and saves them to the class. The trajectory should be aligned previously. You can do this by using the wrapper function that does both steps in one function which is preffered way of doing extraction of waters.

        Parameters
        ----------
        traj : str, optional
            name of the file containing aligned trajectory, by default "aligned_trajectory.xtc"
        topology : str, optional
            topology file, for gromacs give tpr file if you want center of mass centering (important for non constrained runs), by default None
        dist : float, optional
            distance from center of active site for extracting relvenat waters, by default 12.0
        SOL : str, optional
            residue name of water, by default "SOL"
        OW : str, optional
            name of oxygen water, by default "OW"
        HW1 : str, optional
            name of 1st water hydrogen, by default "HW1"
        HW2 : str, optional
            name of 2nd water hydrogen, by default "HW2"

        See Also
        ----------
        select_relevant_waters : use this function instead usually
        """
        # think about another method of selecting waters - ie within some distance of an atom index?
        if topology:
            u = mda.Universe(topology, traj)
        else:
            u = mda.Universe(traj)
        self.nsnaps = len(u.trajectory)
        coordsH = []
        coordsO = []
        # loop over
        for nn, k in enumerate(u.trajectory):
            s = u.select_atoms(self.activesite_aminoacids_mda)
            cc = s.center(None)
            Os = u.select_atoms(
                "name "
                + str(OW)
                + " and resname "
                + str(SOL)
                + " and point "
                + str(cc[0])
                + " "
                + str(cc[1])
                + " "
                + str(cc[2])
                + " "
                + str(dist)
            )
            for i in Os.indices:
                Hs = u.select_atoms(
                    "(name "
                    + str(HW1)
                    + " or name "
                    + str(HW2)
                    + ") and around 1.0 index "
                    + str(i)
                )
                if len(Hs) != 2:
                    print(k, i, len(Hs), "CRAZY WATER WITH WRONG NUMBER OF HYDROGENS")
                for j in Hs.positions:
                    coordsH.append(j)
            # implement check for this if ther is not 2 hhydrogens, also think how to fix around 1.0 to be bonded instead!
            for i, j in zip(Os.positions, Os.indices):
                coordsO.append(i)
        self.X = np.asarray(coordsO)
        self.init_hydrogen(np.asarray(coordsH))
        self.save_H2O()

    def align_mda(
        self,
        align_selection: str = "active_site",
        unaligned_trj_file: str = "unaligned_trajectory.xtc",
        topology: str = None,
        pdb_to_align_to: str = "aligned.pdb",
        output_trj_file: str = "aligned_trajectory.xtc",
    ):
        """
        Align the given trajectory with given snapshot in pdb form to produce aligned trajectory. It is preffered to use the function align_trajectory which generates all necessery files. Uses MDAnalyisis align to align. Probis often produces slightly better alignment.

        Parameters
        ----------
        align_selection : str, optional
            selection to align with. Possible choices are "active_site" or "protein" or a mda analysis selection string. "active_site" uses aminoacids of active site, "protein" uses the whole protein and any other string uses mda selection for alignment selection, by default "active_site"
        unaligned_trj_file : str, optional
            name of the file containing unaligned trajectory, by default "unaligned_trajectory.xtc"
        topology : str, optional
            name of the topology file. Use .tpr for gromacs if you want center of mass centering - very recomended, by default None
        pdb_to_align_to : str, optional
            pdb file name to which trajectory will be aligned to, by default "aligned.pdb"
        output_trj_file : str, optional
            name of the output aligned trajectory file, by default "aligned_trajectory.xtc"

        See Also
        ----------
        align_probis : alignment using probis algorithm
        align_trajectory : use this function instead probably
        """
        if align_selection == "active_site":
            sele = self.activesite_aminoacids_mda
        elif align_selection == "protein":
            sele = "protein"
        else:
            sele = align_selection
        if topology:
            mob = mda.Universe(topology, unaligned_trj_file)
        else:
            mob = mda.Universe(unaligned_trj_file)
        ref = mda.Universe(pdb_to_align_to)
        mda.analysis.align.AlignTraj(
            mob, ref, select=sele, filename=output_trj_file
        ).run()
        del mob
        del ref

    def align_probis(
        self,
        align_selection: str = "active_site",
        probis_exec: str = None,
        unaligned_trj_file: str = "unalignd_trajectory.xtc",
        pdb_to_align_to: str = "aligned.pdb",
        output_trj_file: str = "aligned_trajectory.xtc",
        topology: str = None,
    ):
        """
        Align the given trajectory with given snapshot in pdb form to produce aligned trajectory. It is preffered to use the function align_trajectory which generates all necessery files. Uses probis align to align. See insilab.org/probis-algorithm/ for more info. Probis often produces slightly better alignment then MDAnalyisis.

        Parameters
        ----------
        align_selection : str, optional
            selection to align with. Possible choices are "active_site" or "protein" or a mda analysis selection string. "active_site" uses aminoacids of active site, "protein" uses the whole protein and any other string uses probis selection for alignment selection, by default "active_site"
        probis_exec : str, optional
            relative path of probis executable, if None downloads linux/mac executable from net and runs it in current dir, by default None
        unaligned_trj_file : str, optional
            name of the file containing unaligned trajectory, by default "unaligned_trajectory.xtc"
        topology : str, optional
            name of the topology file. Use tpr for gromacs if you want to center the simulation box to center of mass (important for non constrained runs), by default None
        pdb_to_align_to : str, optional
            pdb file name to which trajectory will be aligned to, by default "aligned.pdb"
        output_trj_file : str, optional
            name of the output aligned trajectory file, by default "aligned_trajectory.xtc"

        """
        if probis_exec == None:
            pwd = os.path.abspath(os.getcwd())
            probis_web = "http://insilab.org/files/probis-algorithm/probis"
            if not (os.path.isfile(pwd + "/probis")):
                import wget

                wget.download(probis_web, pwd + "/probis")
            probis_exec = "probis"
            st = os.stat(probis_exec)
            os.chmod(probis_exec, st.st_mode | stat.S_IEXEC)
        if align_selection == "active_site":
            sele = self.activesite_aminoacids.replace("+", ",")
            sele = (
                ' -motif1 "[:A and (' + sele + ')]" -motif2 "[:A and (' + sele + ')]" '
            )
        elif align_selection == "protein":
            sele = ""
        else:
            sele = align_selection
        if topology:
            u = mda.Universe(topology, unaligned_trj_file)
        else:
            u = mda.Universe(unaligned_trj_file)
        mobile_pdb = "mobile.pdb"
        u.atoms.segments.segids = "A"
        forwrite = u.atoms
        with mda.Writer(output_trj_file, multiframe=True) as W:
            for snap in u.trajectory:
                with mda.Writer(mobile_pdb) as mw:
                    mw.write(forwrite)
                # string to run
                runstring = (
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
                os.system(runstring)
                # read in rota.pdb files and remove all other files
                newal = mda.Universe(
                    pdb_to_align_to[:4] + "A_" + mobile_pdb[:4] + "A.0.rota.pdb"
                )
                # remove aligned pdb files
                removestring = "rm *.rota.pdb probis*.tmp " + mobile_pdb
                os.system(removestring)
                # add cell dimensions and other missing info?
                newal.dimensions = snap.dimensions
                # append snapshot to new trajectory
                selnew = newal.atoms
                W.write(selnew)

    def align_trajectory(
        self,
        traj_fname: str,
        top_fname: str = None,
        every: int = 1,
        align_mode: str = "probis",
        align_target: int = None,
        align_selection: str = "active_site",
        probis_exec: str = None,
        unaligned_trj_file: str = "unaligned_trajectory.xtc",
        output_trj_file: str = "aligned_trajectory.xtc",
    ):
        """
        Function that aligns trajectory with selected method from given trajectory file for given selection. The prefered use is to use the select_relevant_waters function instead of directly calling this.

        Parameters
        ----------
        traj_fname : str
            filename of trajectory.
        top_fname : str, optional
            topology filename. For gromacs use tpr if you want centering to center of mass which is important for non constrained runs. If None no topology is used (this works for pdb trajectories), by default None
        every : int, optional
            select every [every] snapshot from read trajectory, by default 1 - select every snapshot for further analysis
        align_mode : str, optional
            align algorithm to use - either "probis" or "mda" for MDAnalysis, by default "probis"
        align_target : optional
            selects the target snapshot for alignment. Default is the last snapshot of trajectory. If int is given then that frame is taken from FULL trajectory. If string is given then alignment configuration is read from that file (prefer pdb), by default None - last snapshot
        align_selection : str, optional
            selection for alignment. Options: "active_site" for aminoacids in active site, "protein" for whole protein and any other string that can be interpreted by MDAnalyisis or probis for selection depending on align_mode, by default "active_site"
        probis_exec : str, optional
            relative path of probis executable, if None downloads linux/mac executable from net and runs it in current dir. Only needed if align_mode is "probis", by default None
        unaligned_trj_file : str, optional
            name of the file containing unaligned trajectory, by default "unaligned_trajectory.xtc"
        output_trj_file : str, optional
            name of the output aligned trajectory file, by default "aligned_trajectory.xtc"

        See Also
        ---------
        select_relevant_waters
        """
        import MDAnalysis.transformations as trans

        if top_fname:
            mob = mda.Universe(top_fname, traj_fname)
            ref = mda.Universe(top_fname, traj_fname)
        else:
            mob = mda.Universe(traj_fname)
            ref = mda.Universe(traj_fname)
        # center the box to center of mass
        if top_fname.upper().endswith(
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
            )
        ) or traj_fname.upper().endswith(
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
            )
        ):
            protein1 = ref.select_atoms("protein")
            not_protein1 = ref.select_atoms("not protein")
            protein2 = mob.select_atoms("protein")
            not_protein2 = mob.select_atoms("not protein")
            transforms1 = [
                trans.unwrap(protein1),
                trans.center_in_box(protein1, wrap=True),
                trans.wrap(not_protein1),
            ]
            transforms2 = [
                trans.unwrap(protein2),
                trans.center_in_box(protein2, wrap=True),
                trans.wrap(not_protein2),
            ]
            ref.trajectory.add_transformations(*transforms1)
            mob.trajectory.add_transformations(*transforms2)
        else:
            print(
                "cannot center the protein - consider using tpr, prmtop, pdb or other formats that contain bond information\n Continuing..."
            )

        self.nsnaps = len(mob.trajectory)
        ref.atoms.segments.segids = "A"
        if align_target:
            if type(align_target) == int:
                wr = ref.atoms
                wr.write("aligned.pdb", frames=ref.trajectory[[align_target]])
                align_target = "aligned.pdb"
        else:
            wr = ref.atoms
            wr.write("aligned.pdb", frames=ref.trajectory[[-1]])
            align_target = "aligned.pdb"
        self.align_file = align_target
        if every > 1:
            # for writing to pdb file perhaps make an wrapper that supresses stupid warnings
            # silence warnings with:
            # import warnings
            # warnings.filterwarnings("ignore", message="Numerical issues were encountered ")
            # #### or ####
            # with warnings.catch_warnings():
            #   warnings.simplefilter('ignore')
            with mda.Writer(unaligned_trj_file, multiframe=True) as W:
                for i in mob.trajectory[::every]:
                    W.write(mob.atoms)
        else:
            shutil.copy2(traj_fname, unaligned_trj_file)
        if align_mode == "probis":
            self.align_probis(
                align_selection=align_selection,
                topology=top_fname,
                probis_exec=probis_exec,
                unaligned_trj_file=unaligned_trj_file,
                pdb_to_align_to=align_target,
                output_trj_file=output_trj_file,
            )
        elif align_mode == "mda":
            self.align_mda(
                align_selection=align_selection,
                output_trj_file=output_trj_file,
                topology=top_fname,
                unaligned_trj_file=unaligned_trj_file,
                pdb_to_align_to=align_target,
            )

    def select_relevant_waters(
        self,
        traj_fname: str,
        top_fname: str = None,
        every: int = 1,
        dist: float = 12.0,
        align_mode: str = "probis",
        align_target=None,
        align_selection: str = "active_site",
        probis_exec: str = None,
        SOL: str = "SOL",
        OW: str = "OW",
        HW1: str = "HW1",
        HW2: str = "HW2",
    ):
        """
        Function for selecting relevant waters from given trajectory and selection. The trajectory is first aligned and waters are selected from aligned trajectory. If you do not wish to align the trajectory use extract_waters_from_trajectory.

        Parameters
        ----------
        traj_fname : str
            file name containing trajectory
        top_fname : str, optional
            topology filename. For gromacs use tpr if you want centering to center of mass which is important for non constrained runs. If None no topology is used (this works for pdb trajectories), by default None
        every : int, optional
            select every [every] snapshot from read trajectory, by default 1 - select every snapshot for further analysis
        dist : float, optional
            distance around centre of active site inside which waters will be searched for, by default 12.0
        align_mode : str, optional
            align algorithm to use - either "probis" or "mda" for MDAnalysis, by default "probis"
        align_target : optional
            selects the target snapshot for alignment. Default is the last snapshot of trajectory. If int is given then that frame is taken from FULL trajectory. If string is given then alignment configuration is read from that file (prefer pdb), by default None - last snapshot
        align_selection : str, optional
            selection for alignment. Options: "active_site" for aminoacids in active site, "protein" for whole protein and any other string that can be interpreted by MDAnalyisis or probis for selection depending on align_mode, by default "active_site"
        probis_exec : str, optional
            relative path of probis executable, if None downloads linux/mac executable from net and runs it in current dir. Only needed if align_mode is "probis", by default None
            , by default None
        SOL : str, optional
            residue name of water, by default "SOL"
        OW : str, optional
            name of oxygen water, by default "OW"
        HW1 : str, optional
            name of 1st water hydrogen, by default "HW1"
        HW2 : str, optional
            name of 2nd water hydrogen, by default "HW2"

        """
        self.align_trajectory(
            traj_fname=traj_fname,
            top_fname=top_fname,
            every=every,
            align_mode=align_mode,
            align_target=align_target,
        )
        self.extract_waters_from_trajectory(
            dist=dist, topology=top_fname, SOL=SOL, OW=OW, HW1=HW1, HW2=HW2
        )

    def init_hydrogen(self, Hcoords: np.ndarray):
        """
        Saves hydrogen coordinates to class array for further analysis and produces orientation vectors.

        Parameters
        ----------
        Hcoords : np.ndarray
            numpy array containing xyz coordinates of hydrogens. Two hydrogens correspoding to same oxygen should be successive in this array.

        """
        self.H = Hcoords
        self.H1 = self.H[::2, :]
        self.H2 = self.H[1:, :]
        self.H2 = self.H2[::2, :]
        self.v1 = []
        self.v2 = []
        for o, h1, h2 in zip(self.X, self.H1, self.H2):
            a = h1 - o
            b = h2 - o
            if (np.linalg.norm(h1 - o) > 1.2) or (np.linalg.norm(h2 - o) > 1.2):
                print("bad input")
                sys.exit()
            self.v1.append(a)
            self.v2.append(b)
        self.v1 = np.asarray(self.v1)
        self.v2 = np.asarray(self.v2)

    def cluster(
        self,
        maxsamp: int = None,
        clustering_type: str = "type_by_type",
        options: list = None,
        clustering_algorithm: str = "OPTICS",
        verbose: int = 0,
        debugO: int = 0,
        debugH: int = 0,
        plot_reachability: int = 0,
        plotend: int = 0,
    ):
        """
        Main clustering function. Invoke to perform clustering analysis of water molecules within previously selected region.

        Parameters
        ----------
        maxsamp : int, Optional, default: None - defaults to number of snapshots.
            maxsamp argument for clustering analysis. Maximum number to be given should be number of snapshots used to produce water molecules data. In practice set it to maximum if type_by_type is used.
        clustering_type: string, Optional, default:"type_by_type"
            Selects clustering strategy. "type_by_type" works best. Possible choices:
                type_by_type - First checks for conserved, then half conserved and finaly semi conserved waters using hydrogen orientation analysis. See scan_oxygen_clustering for more info.
                converge_oxygen_clustering - Tries to converge clustering of oxygens before trying to identify water types using hydrogen orientation analysis.
                loop_conserved_then_converge - finds conserved waters using scan_oxygen_clustering methodology and then uses cinverge_oxygen_clustering to find the rest of the waters.
                oxygen_clustering - uses scan_oxygen_clustering to find one type of conserved water. See scan_oxygen_clustering.
                cluster_once - performes clustering once with given minsap for options
        options: various, Optional, default: None
            Aditional options for different clustering types, if it is None default option is used:
                type_by_type : lower_minsamp - float, lowest minsamp at which scanning is stopped as percentage of maxsamp
                loop_conserved_then_converge and converge_oxygen_clustering : list of two floats - [xi, treshold for convergence]
                oxygen_clustering : type of water can be a list containing any combo of "FCW", "HCW" or "WCW" - str; if "onlyO" only oxygen clustering is reported
                cluster_once: list containing minsamp, xi and list which contains any combination of "FCW", "HCW" or "WCW"; or "onlyO"
        clustering_algorithm: str, Optional, default: "OPTICS"
            Type of clustering algorithm for oxygen clustering. OPTICS seems to give best results. Other choices: DBSCAN and HDBSCAN
        verbose: int, Optional, default : 0
            0 only print main loop envts; 1 - print for selected waters; 2 - all considered
        debugO: int, Optional, default : 0
            for oxygen clustering:0 - plot nothing; 1 - plot only selected; 2 - plot all considered
        debugH: int, Optional, default : 0
            for hydrogen orient :0 - plot nothing; 1 - plot only selected; 2 - plot all considered
        plotter: int, Optional, default : 0
            0 - plot only clustering; 1 plot reachability as well where applicable
        plotend: int, Optional, default : 0
            0 plots when found/stopped; 1 - plots when all are found; only put 1 if both debugH and debugO are <2

        """
        if len(self.X) == 0 or len(self.v1) == 0 or len(self.v2) == 0:
            print("oxygen and hydrogen arrays not defined - STOPPED")
            return
        if self.target_oxygen_cluster_size is None:
            if self.nsnaps is None:
                if maxsamp is None:
                    print("maxsamp not defined - STOPPED")
                    return
                else:
                    self.maxsamp = maxsamp
                    self.target_oxygen_cluster_size = self.maxsamp
            else:
                self.target_oxygen_cluster_size = self.nsnaps
        if maxsamp == None:
            if self.nsnaps:
                maxsamp = self.nsnaps
            else:
                print("no maxsamp given - STOPPED")
                return
        self.maxsamp = maxsamp
        self.verbose = verbose
        self.debugO = debugO
        self.debugH = debugH
        self.plotter = plot_reachability
        self.plotend = plotend
        self.clustering_type = clustering_type
        self.clustering_options = options
        self.clustering_algorithm = clustering_algorithm
        if clustering_algorithm == "OPTICS":
            self.isHDBSCAN = False
            self.isDBSCAN = False
        if clustering_algorithm == "DBSCAN":
            self.isHDBSCAN = False
            self.isDBSCAN = True
        if clustering_algorithm == "HDBSCAN":
            self.isHDBSCAN = True
            self.isDBSCAN = False
        self.save_clustering_options()
        if clustering_type == "cluster_once":
            self.delete_after_find = False
            if options:
                self.recalculate_reachability(options[0])
                self.oxygen_clustering(
                    minsamp=options[0], xi=options[1], whichH=options[2]
                )
            else:
                self.recalculate_reachability(int(self.target_oxygen_cluster_size / 8))
                self.oxygen_clustering(
                    minsamp=int(self.target_oxygen_cluster_size / 8), xi=0.01
                )
        if clustering_type == "type_by_type":
            self.loop_type_by_type(lower_minsamp=options)
        if clustering_type == "oxygen_clustering":
            if options:
                self.scan_oxygen_clustering(whichH=[options])
            else:
                self.scan_oxygen_clustering()
        if clustering_type == "converge_oxygen_clustering":
            if options:
                self.converge_oxygen_clustering(xi=options[0], trsh=options[1])
            else:
                self.converge_oxygen_clustering()
        if clustering_type == "loop_conserved_then_converge":
            if options:
                self.loop_conserved_then_converge_rest(xi=options[0], trsh=options[1])
            else:
                self.loop_conserved_then_converge_rest()
        self.save_clustering_options()
        self.save_results()
        # restore defaults
        self.isDBSCAN = False
        self.isHDBSCAN = False
        self.verbose = 1
        self.debugO = 0
        self.debugH = 0
        self.plotter = 0
        self.plotend = 0
        self.maxsamp = None
        self.clustering_type = None
        self.clustering_options = None
        self.clustering_algorithm = None

    def save_clustering_options(self, fname: str = "clust_options.dat"):
        """
        Function for saving clustering options.
        Saves clustering options for restart.

        Parameters
        ----------
        fname : str, optional
            filename where clustering options shall be saved, by default "clust_options.dat"
        """
        """

        """
        with open(fname, "w") as f:
            print(self.maxsamp, file=f)
            print(self.clustering_algorithm, file=f)
            print(self.clustering_type, file=f)
            print(self.verbose, file=f)
            print(self.debugO, file=f)
            print(self.debugH, file=f)
            print(self.plotter, file=f)
            print(self.plotend, file=f)
            if isinstance(self.clustering_options, Iterable):
                for i in self.clustering_options:
                    print(i, file=f)
            else:
                print(self.clustering_options, file=f)

    def restart_cluster(
        self,
        options_file="clust_options.dat",
        water_file="water_coords_restart.dat",
        results_file="Clustering_results_temp.dat",
        run=True,
    ):
        """
        Function which restarts clustering from restart files.

        Replace this function and doc string for your own project

        Parameters
        ----------
        options_file : str, Optional, default: "clust_options.dat"
            Restart file that contains clustering parameters. You shouldn't have to change this.
        water_file : str, Optional, default: "water_coords_restart.dat"
            file containing water coordinates after each deletion.
        resluts_file : str, Optional, default: "Clustering_results_temp.dat"
            file containing water clusters found so far.
        run : bool, default: True
            restarts the clustering if True after loading the clustering options and results.

        """
        if os.path.isfile(water_file):
            self.load_H2O(fname=water_file)
        elif len(self.X) == 0 or len(self.v1) == 0:
            print("No input water coordinates given STOPPING")
            return
        if os.path.isfile(results_file):
            self.read_results(results_file)
        if os.path.isfile(options_file):
            f = open(options_file)
            lines = f.read().splitlines()
            self.maxsamp = int(lines[0])
            self.nsnaps = self.maxsamp + 1
            self.clustering_algorithm = lines[1].strip(" ")
            self.clustering_type = lines[2].strip(" ")
            self.verbose = int(lines[3])
            self.debugO = int(lines[4])
            self.debugH = int(lines[5])
            self.plotter = int(lines[6])
            self.plotend = int(lines[7])
            self.clustering_options = None
            if self.clustering_type == "type_by_type":
                if lines[8].strip(" ") != "None":
                    self.clustering_options = float(lines[8])
            if self.clustering_type == "oxygen_clustering":
                if lines[8].strip(" ") != "None":
                    self.clustering_options = lines[8].strip(" ")
            if self.clustering_type == "converge_oxygen_clustering":
                if lines[8].strip(" ") != "None":
                    self.clustering_options = [float(lines[8]), float(lines[9])]
            if self.clustering_type == "loop_conserved_then_converge":
                if lines[8].strip(" ") != "None":
                    self.clustering_options = [float(lines[8]), float(lines[9])]
            if run:
                self.cluster(
                    maxsamp=self.maxsamp,
                    clustering_type=self.clustering_type,
                    options=self.clustering_options,
                    clustering_algorithm=self.clustering_algorithm,
                    verbose=self.verbose,
                    debugO=self.debugO,
                    debugH=self.debugH,
                    plot_reachability=self.plotter,
                    plotend=self.plotend,
                )

    def save_H2O(self, fname="water_coordinates.dat"):
        """
        saves oxygen xyz coordiantes and hydrogens orientations in 9 columns.

        Parameters
        ----------
        fname : str, optional
            filename to save to, by default "water_coordinates.dat"
        """
        np.savetxt(fname, np.c_[self.X, self.v1, self.v2])

    def load_H2O(self, fname="water_coordinates.dat"):
        """
        loads water coordinates from a file. The water coordinates file should contain 9 columns: xyz coordinates of oxygen, and xyz orientations of hydrogen 1 and 2 for given oxygen

        Parameters
        ----------
        fname : str, optional
            filename to read from, by default "water_coordinates.dat"

        """
        data = np.loadtxt(fname)
        self.X = data[:, :3]
        self.v1 = data[:, 3:6]
        self.v2 = data[:, 6:9]

    def scan_oxygen_clustering(self, whichH: List[str] = ["FCW", "HCW", "WCW"]):
        """
        Main loop - loops over oxygen clustering parameter space (minsamp and xi) and calculate reachability - if a clustering with satisfactory oxygen clustering and hydrogen orientation clustering (optional) is found elements of that oxygen cluster are removed from the data set and oxygen clustering starts from the beggining. loops until no satisfactory clusterings are found.
        lowest minsamp is defined by self.lower_minsamp

        Parameters
        ----------
        whichH : List[str], optional
            which water types to scan for (FCW, HCW or WCW), default - all

        """
        self.setdebuggers()
        found = True
        print("Searching for waters")
        # set of xi's that will be used for oxygen OPTICS clustering
        xis = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00001]
        if len(self.X) < self.maxsamp:
            found = False
        while found:
            found = False
            # loop over minsamps- from N(snapshots) to 0.75*N(snapshots)
            for i in reversed(
                range(int(self.maxsamp * self.lower_minsamp), self.maxsamp + 1)
            ):
                print("minsamp", i)
                # recalculate reachability - OPTICS reachability has to be recaculated when changing minsamp
                self.recalculate_reachability(i)
                # loop over xi
                if self.isHDBSCAN:
                    j = 0
                    waters = self.oxygen_clustering(i, j, whichH)
                    if len(waters) > 0:
                        found = True
                        break
                if self.isDBSCAN and not (self.isHDBSCAN):
                    xis = np.linspace(
                        np.min(self.reachability),
                        np.mean(self.reachability[np.isfinite(self.reachability)]),
                        21,
                    )[1:]
                if not (self.isHDBSCAN):
                    for j in xis:
                        waters = self.oxygen_clustering(i, j, whichH)
                        if len(waters) > 0:
                            found = True
                            break
                if found:
                    break
            # check if size of remaining data set is bigger then maxsamp
            if len(self.X) < self.maxsamp:
                found = False
        print("done")
        if (self.debugH == 1 or self.debugO == 1) and self.plotend == 1:
            plt.show()

    def converge_oxygen_clustering(self, xi: float = 0.001, trsh: int = None):
        """
        Check convergence of oxygen clustering wrt to minsamp for given xi. After clustering converges, calculate conserved water molecules. Convergence is met if two successive oxygen clustering have same number of clusters, and sum of differences of cluster sizes is smaller then given treshold

        Parameters
        ----------
        xi : float, optional
            OPTICS parameter xi or eta for DBSCAN, by default 0.001
        trsh : int, optional
            convergence treshold. Number represents the sum of differences between different cluster sizes as compared to previous iteration of clustering preformed with smaller minsamp, by default None - clusters need to be equal size
        """
        self.setdebuggers()
        print("finding converging oxygen clustering")
        bestdiff = len(self.X)
        bestminsamp = self.maxsamp
        target_minsamp = 0
        print("treshold", trsh, "xi", xi)
        self.delete_after_find = False
        # loop over minsamps- from 0.25*N(snapshots) to N(snapshots)
        for i, minsamp in enumerate(range(int(self.maxsamp * 0.25), self.maxsamp + 1)):
            print("minsamp", minsamp)
            # recalculate reachability - OPTICS reachability has to be recaculated when changing minsamp
            self.recalculate_reachability(minsamp)
            # calculate optics clustering from reachability with wanted xi that has been previously calculated
            if self.isDBSCAN:
                labels = cluster_optics_dbscan(
                    reachability=self.clust.reachability_,
                    ordering=self.clust.ordering_,
                    core_distances=self.clust.core_distances_,
                    eps=xi,
                )
            else:
                labels = cluster_optics_xi(
                    reachability=self.clust.reachability_,
                    predecessor=self.clust.predecessor_,
                    ordering=self.clust.ordering_,
                    min_samples=minsamp,
                    xi=xi,
                )[0]
            # Calculate number of elements in each cluster and sort them according to size
            (values, counts) = np.unique(labels, return_counts=True)
            counts = counts[np.argsort(counts)]
            if i == 0:
                oldlabels = labels
                oldvalues = values
                oldcounts = counts
            else:
                converged = False
                # check number of clusters
                if len(values) == len(oldvalues):
                    converged = True
                # check if cluster sizes match
                if converged:
                    print("samesize")
                    diff = 0
                    for i, j in zip(counts, oldcounts):
                        diff += abs(i - j)
                    print("diff:", diff)
                    if trsh:
                        print("TRSH", trsh)
                        if diff > trsh:
                            converged = False
                    else:
                        if diff < bestdiff:
                            bestdiff = diff
                            bestminsamp = minsamp
                        converged = False
                oldlabels = labels
                oldcounts = counts
                oldvalues = values
                # if converged label ideal minsamp
                if converged:
                    print("convergence criteria met")
                    target_minsamp = minsamp
                    break
        if not (converged):
            print("best minsamp diff", bestminsamp)
            target_minsamp = bestminsamp
        print("done")
        if target_minsamp > 0 or not (trsh):
            waters = self.oxygen_clustering(target_minsamp, xi)
            if len(waters) > 0:
                self.add_visualisation_waters(waters)
            if (self.debugH == 1 or self.debugO == 1) and self.plotend == 1:
                plt.show()

    def loop_conserved_then_converge_rest(self, xi: float = 0.001, trsh: int = None):
        """
        Combines scanning oxygen clustering parameters for conserved waters until all conserved waters are found. Then tries converging the clustering for remaining data and calculate half and semi conserved.

        Parameters
        ----------
        xi : float, optional
            OPTICS parameter xi or eta for DBSCAN, by default 0.001
        trsh : int, optional
            convergence treshold. Number represents the sum of differences between different cluster sizes as compared to previous iteration of clustering preformed with smaller minsamp, by default None - clusters need to be equal size

        """
        if self.isHDBSCAN:
            self.plotter = 0
        self.lower_minsamp = 0.5
        self.scan_oxygen_clustering(whichH=["FCW"])
        self.converge_oxygen_clustering(xi=xi, trsh=trsh)

    def loop_type_by_type(self, lower_minsamp: float = 0.25):
        """
        Scans oxygen clustering parameters for conserved waters until all conserved waters are found and deletes found entries from the database. Then does the same for half conserved, and then semi conserved respectively.

        Parameters
        ----------
        lower_minsamp : float, optional
            lowest percentage of total number of snapshots to be used for minsamp parameter of clustering, by default 0.25

        """
        if self.isHDBSCAN:
            self.plotter = 0
        self.lower_minsamp = 0.25 if lower_minsamp is None else lower_minsamp
        self.scan_oxygen_clustering(whichH=["FCW"])
        self.scan_oxygen_clustering(whichH=["HCW"])
        self.scan_oxygen_clustering(whichH=["WCW"])

    def calc_clustering(self, minsamp=None, xi=0.001, whichH=["FCW", "HCW", "WCW"]):
        """Calculates waters using OPTICS clustering for given minsamp,xi and hydrogen types
        Marked for depraction? Could be used for testing?
        """
        if not (minsamp):
            minsamp = self.maxsamp * 0.8
        # recalculate reachability - OPTICS reachability has to be recaculated when changing minsamp
        self.recalculate_reachability(minsamp)
        # calculate optics clustering from reachability with wanted xi that has been previously calculated
        labels = cluster_optics_xi(
            reachability=self.clust.reachability_,
            predecessor=self.clust.predecessor_,
            ordering=self.clust.ordering_,
            min_samples=minsamp,
            xi=xi,
        )[0]
        waters = self.oxygen_clustering(minsamp, xi, whichH)
        if len(waters) > 0:
            self.add_visualisation_waters(waters)
        if (self.debugH == 1 or self.debugO == 1) and self.plotend == 1:
            plt.show()

    def oxygen_clustering(
        self, minsamp: int, xi: float, whichH: str = ["FCW", "HCW", "WCW"]
    ):
        """
        Analyses oxygen OPTICS clustering based on minsamp and xi. For oxygen clusters which have the size around number of samples, the hydrogen orientation analysis is performed and type of water molecule and coordinates are given back.

        Parameters
        ----------
        minsamp: int
            minimum number of samples for OPTICS clustering
        xi: float
            xi for OPTICS clustering, eta for DBSCAN
        whichH: list of str
            list of which hydrogen types to analyse (FCW, HCW or WCW) or only O clustering (onlyO) - NOT IN list

        Returns
        -------
        waters : List
            returns list which contains a list of positions for O,H,H atom coordinates and type of water
        """
        waters = []
        # calculate optics clustering from reachability with wanted xi that has been previously calculated
        if self.isDBSCAN:
            clusters = cluster_optics_dbscan(
                reachability=self.clust.reachability_,
                ordering=self.clust.ordering_,
                core_distances=self.clust.core_distances_,
                eps=xi,
            )
        elif self.isHDBSCAN:
            clusters = self.clust.labels_
        else:
            clusters = cluster_optics_xi(
                reachability=self.clust.reachability_,
                predecessor=self.clust.predecessor_,
                ordering=self.clust.ordering_,
                min_samples=minsamp,
                xi=xi,
            )[0]
        # Debug stuff
        if self.debugO > 0:
            ff = self.debug_O(minsamp, xi, clusters)
        # Loop over all oxygen clusters (-1 is non cluster)
        for k in np.sort(np.unique(clusters[clusters != -1])):
            # Number of elements in oxygen cluster
            neioc = np.count_nonzero(clusters == k)
            # If number of elements in oxygen cluster is  Nsnap*0.85<Nelem<Nsnap*1.15 then ignore
            if (
                neioc < self.target_oxygen_cluster_size * (2 - self.numbpct)
                and neioc > self.target_oxygen_cluster_size * self.numbpct
            ):
                if self.verbose > 0:
                    print(f"O clust {k}, size {len(clusters[clusters==k])}\n")
                O = np.mean(self.X[clusters == k], axis=0)
                water = [O]
                if not (whichH[0] == "onlyO"):
                    # Construct array of hydrogen orientations
                    orientations = []
                    for i in np.argwhere(clusters == k):
                        orientations.append(self.v1[i])
                        orientations.append(self.v2[i])
                    orientations = np.asarray(orientations)[:, 0, :]
                    # Analyse clustering with hydrogen orientation analysis and more debug stuff
                    hyd = hydrogen_orientation_analysis(
                        orientations,
                        self.numbpct_hyd_orient_analysis,
                        self.kmeansangcutoff,
                        self.kmeansinertiacutoff,
                        self.conserved_optics_angdiff_cutoff,
                        self.conserved_optics_angstd_cutoff,
                        self.other_waters_minsamp_pct,
                        self.half_conserved_angle_limit,
                        self.conserved_optics_angstd_cutoff,
                        self.wiggly_conserved_std_limit,
                        self.wiggly_explained,
                        njobs=self.njobs,
                        verbose=self.verbose,
                        debugH=self.debugH,
                        plotreach=(self.plotter == 1),
                        which=whichH,
                    )
                    if self.plotend == 0 and self.debugH > 0:
                        plt.show()
                    if len(hyd) > 0:
                        # add water atoms for pymol visualisation
                        for i in hyd:
                            water = [O]
                            water.append(O + i[0])
                            water.append(O + i[1])
                            water.append(i[2])
                            waters.append(water)
                        # debug
                        if self.debugO == 1 and self.plotend == 0 and self.debugH == 0:
                            plt.show()
                        if self.delete_after_find:
                            self.delete_data(np.argwhere(clusters == k))
                            break
                    elif self.debugO == 1:
                        plt.close(ff)
                else:
                    print("O_clust")
                    water.append(O)
                    water.append(O)
                    water.append("O_clust")
                    waters.append(water)
                    if self.delete_after_find:
                        self.delete_data(np.argwhere(clusters == k))
                        break
        if self.debugO == 1 and len(waters) == 0:
            plt.close(ff)
        if len(waters) > 0:
            self.add_visualisation_waters(waters)
        return waters

    def delete_data(self, elements):
        """deletes the data for oxygen and hydrogen"""
        self.X = np.delete(self.X, elements, 0)
        self.v1 = np.delete(self.v1, elements, 0)
        self.v2 = np.delete(self.v2, elements, 0)
        self.save_H2O(fname="water_coords_restart1.dat")
        os.rename("water_coords_restart1.dat", "water_coords_restart.dat")
        self.save_results(fname="Clustering_results_temp1.dat")
        os.rename("Clustering_results_temp1.dat", "Clustering_results_temp.dat")
        os.rename(
            "Type_Clustering_results_temp1.dat", "Type_Clustering_results_temp.dat"
        )

    def recalculate_reachability(self, msp):
        """Performs OPTICS or HDBSCAN clustering based on min samp and saves in "global" variables"""
        if self.isHDBSCAN:
            self.clust = hdbscan.HDBSCAN(
                min_cluster_size=int(self.target_oxygen_cluster_size * self.numbpct),
                min_samples=int(msp),
                max_cluster_size=int(
                    self.target_oxygen_cluster_size * (2 - self.numbpct)
                ),
                cluster_selection_method="eom",
                algorithm="best",
                core_dist_n_jobs=self.njobs,
            )
        else:
            self.clust = OPTICS(min_samples=int(msp), n_jobs=self.njobs)
        self.clust.fit(self.X)
        self.space = np.arange(len(self.X))
        if self.isHDBSCAN:
            self.labels = self.clust.labels_
        else:
            self.reachability = self.clust.reachability_[self.clust.ordering_]
            self.labels = self.clust.labels_[self.clust.ordering_]

    def debug_O(self, i, j, labels):
        """For debuging oxygen clustering"""
        if self.verbose > 0:
            (aa, bb) = np.unique(labels, return_counts=True)
            so = (
                f"Oxygen clustering minsamp={i} and xi={j}, {len(np.unique(labels[labels!=-1]))} clusters \n"
                f"Required N(elem) range:{self.maxsamp*self.numbpct:.2f} to {self.maxsamp}; (maxsamp={self.maxsamp} and numbpct={self.numbpct:.2f})\n"
                f"N(elements) for each cluster: {bb}\n"
            )
        if self.debugO > 0:
            fig = plt.figure()
            if self.plotter == 1:
                ax = fig.add_subplot(1, 2, 1, projection="3d")
            else:
                ax = fig.add_subplot(1, 1, 1, projection="3d")
            for k in np.sort(np.unique(labels)):
                neioc = np.count_nonzero(labels == k)
                jaba = self.X[labels == k]
                s = 10
                if k == -1:
                    s = 0.25
                ax.scatter(
                    jaba[:, 0],
                    jaba[:, 1],
                    jaba[:, 2],
                    label=f"{k} ({len(labels[labels==k])})",
                    s=s,
                )
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_zlabel("Z")
            ax.set_title(so)
            ax.legend()
            if self.plotter == 1:
                lblls = labels[self.clust.ordering_]
                ax = fig.add_subplot(1, 2, 2)
                plt.gca().set_prop_cycle(None)
                ax.plot(self.space, self.reachability)
                for clst in np.unique(lblls):
                    if clst == -1:
                        ax.plot(
                            self.space[lblls == clst],
                            self.reachability[lblls == clst],
                            label=f"{clst} ({len(self.space[lblls==clst])})",
                            color="blue",
                        )
                    else:
                        ax.plot(
                            self.space[lblls == clst],
                            self.reachability[lblls == clst],
                            label=f"{clst} ({len(self.space[lblls==clst])})",
                        )
                ax.legend()
        if self.verbose > 0:
            print(so)
        if self.debugO == 2 and self.plotend == 0:
            plt.pause(0.01)
        return fig

    def save_results(self, fname="Clustering_results.dat"):
        """
        Saves clustering results. Once saved the file can be used to load the results or restart the clustering run or generate visualisation.

        Parameters
        ----------
        fname : str, optional
            file name primer for saving the results, by default "Clustering_results.dat"
        """
        if not (self.analyse_hydrogen):
            np.savetxt(fname, np.c_[self.waterO])
        else:
            np.savetxt(fname, np.c_[self.waterO, self.waterH1, self.waterH2])
        np.savetxt("Type_" + fname, np.c_[self.water_type], fmt="%s")

    def read_results(self, fname="Clustering_results.dat"):
        """
        Reads clustering results from files. After the results are read the clustering can be continued(restarted), or one can visualse the results.

        Parameters
        ----------
        fname : str, optional
            file name primer of the results, by default "Clustering_results.dat"
        """
        self.water_type = []
        self.waterO = []
        self.waterH1 = []
        self.waterH2 = []
        coords = np.loadtxt(fname)
        if coords.shape[1] == 3:
            for i in coords:
                self.waterO.append(i)
        else:
            O = coords[:, :3]
            H1 = coords[:, 3:6]
            H2 = coords[:, 6:9]
            for i, j, k in zip(O, H1, H2):
                self.waterO.append(i)
                self.waterH1.append(j)
                self.waterH2.append(k)
        types = np.loadtxt("Type_" + fname, dtype=str)
        for i in types:
            self.water_type.append(i)

    def visualise_pymol_old(self, pdbfile: str):
        """Marked for depraction/ revision to load the outputed pdb file"""
        self.visualise_init(pdbfile, self.activesite_aminoacids)
        self.visualise_waters()

    def visualise_init(self, fname, am):
        import pymol
        from pymol import cmd

        pymol.finish_launching(["pymol", "-q"])
        cmd.load(fname)
        self.origin_water = cmd.get_unused_name("backup_water_")
        cmd.create(self.origin_water, "bymolecule (solvent and index 1)")
        cmd.hide("everything")
        # naredi samo polimer za def povšine
        tmpObj = cmd.get_unused_name("_tmp")
        cmd.create(tmpObj, "( all ) and polymer", zoom=0)
        # amino kiseline u aktivnom mestu
        self.aminokis_u_am = cmd.get_unused_name("ak_u_am_")
        cmd.select(self.aminokis_u_am, " resi " + am + " and polymer ")
        cmd.show("licorice", self.aminokis_u_am)
        # dodaj pseudoatom u centru aktivnog mesta
        imezapseudoatom = cmd.get_unused_name("centaraktivnogmesta_")
        self.active_site_center = imezapseudoatom
        cmd.pseudoatom(imezapseudoatom, self.aminokis_u_am)
        # sredi pseudoatom
        cmd.hide(representation="everything", selection=imezapseudoatom)
        # cmd.show("sphere",imezapseudoatom)
        # cmd.set ("sphere_scale",0.1,imezapseudoatom)
        cmd.center(imezapseudoatom)
        # napravi povrsinu proteina
        protein = cmd.get_unused_name("samo_protein_")
        cmd.select(protein, "polymer")
        povrsina = cmd.get_unused_name("povrsina_protein_")
        cmd.create(povrsina, protein)
        cmd.show("surface", povrsina)
        cmd.color("gray70", povrsina)
        cmd.set("transparency", 0.5, povrsina)
        cmd.center(imezapseudoatom)
        cmd.delete(fname[: fname.rfind(".")])
        cmd.save("clustering_session.pse")

    def draw_waters(self, O, H1, H2, wname):
        import pymol
        from pymol import cmd

        # copy the water
        cmd.copy(wname, self.origin_water)
        # change coordinates of oxygen atom from water
        cmd.alter_state(
            0, wname, "(x,y,z)=(" + str(O[0]) + "," + str(O[1]) + "," + str(O[2]) + ")"
        )
        ##OPTION ONE - rotate to match dipole vectors after rotationthe plane of the molecule has to be rotated around dipoleaxis to match the plane of the cluster calculated water
        # perhaps cmd.rotate or below
        # cmd.alter_state("0",wname+" and elem O", "(x,y,z)=("+str([0])+","+str(O[1])+","+str(O[2])+")")
        indeciesH = []
        cmd.iterate_state(
            -1,
            wname + " and elem H",
            "indeciesH.append(index)",
            space={"indeciesH": indeciesH},
        )
        cmd.alter_state(
            0,
            wname + " and index " + str(indeciesH[0]),
            "(x,y,z)=(" + str(H1[0]) + "," + str(H1[1]) + "," + str(H1[2]) + ")",
        )
        cmd.alter_state(
            0,
            wname + " and index " + str(indeciesH[1]),
            "(x,y,z)=(" + str(H2[0]) + "," + str(H2[1]) + "," + str(H2[2]) + ")",
        )
        if str(wname).startswith("Orient"):
            cmd.show("spheres", wname)
            cmd.set("sphere_scale", 0.1, wname)
        else:
            cmd.show("sticks", wname)

    def visualise_waters(self):
        looper = zip(self.waterO, self.waterH1, self.waterH2, self.water_type)
        for O, H1, H2, tip in looper:
            wname = cmd.get_unused_name(tip + "_")
            self.all_waters.append(wname)
            print(wname)
            self.draw_waters(O, H1, H2, wname)
        cmd.delete(self.origin_water)
        cmd.reset()
        cmd.center(self.active_site_center)
        cmd.distance("polar_contacts", self.aminokis_u_am + " or sol", "sol", mode=2)
        cmd.hide("labels")
