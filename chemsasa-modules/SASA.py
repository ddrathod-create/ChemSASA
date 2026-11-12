import numpy as np
import MDAnalysis as mda
import freesasa
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from res_func_grp import BACKBONE, SIDECHAINS

class SASA:
    def __init__(self, topology, trajectory=None, unique_residues=None,
                 classifier_path=None,
                 select="protein and not name H*", output_dir=None):

        self.topology = topology
        self.trajectory = trajectory
        self.has_trajectory = trajectory is not None
        self.classifier_path = classifier_path
        self.select = select
        self.output_dir = Path(output_dir) if output_dir else None

        self.universe = mda.Universe(topology, trajectory) if trajectory else mda.Universe(topology)
        self.atomgroup = self.universe.select_atoms(select)
        self.unique_residues = unique_residues

        self.classifier = freesasa.Classifier(classifier_path)
        self.parameters = freesasa.Parameters({
            'algorithm': freesasa.LeeRichards,
            'probe-radius': 1.4,
            'n-slices': 100
        })

        self.radii = self._compute_radii()

        self.func_grp_labels = sorted(self._functional_groups())
        self._func_grp_index = {g: i for i, g in enumerate(self.func_grp_labels)}

        self.residue_labels = self._build_residue_labels()
        self._residue_index = {key: i for i, key in enumerate(self.residue_labels)}

        self.atom_labels = [(a.segid, a.resid, a.resname, a.name) for a in self.atomgroup]

        # results
        self.total_area = None
        self.residue_area = None
        self.func_grp_area = None
        self.atom_area = None
        self.structure_total_area = None
        self.structure_atom_area = None
        self.structure_residue_area = None
        self.structure_func_grp_area = None


    # ---------------- setup ----------------

    def _compute_radii(self):
        radii, n_unknown = [], 0
        for atom in self.atomgroup:
            r = self.classifier.radius(atom.resname, atom.name)
            if r < 0:
                n_unknown += 1
                r = 1.4
            radii.append(r)
        if n_unknown:
            print(f"Warning: {n_unknown}/{len(radii)} atoms had no classifier match, used r = 1.4 Å")
        return np.array(radii, dtype=np.float64)

    def _functional_groups(self):
        all_groups = set(BACKBONE.values())
        for sidechain in SIDECHAINS.values():
            all_groups.update(sidechain.values())
        return all_groups

    def _build_residue_labels(self):
        return [(res.segid, res.resid, res.resname) for res in self.atomgroup.residues]

    # ---------------- static structure (PDB, single frame only) ----------------

    def calculate_sasa_structure(self):
        """Uses freesasa's own file parser directly on self.topology
        returns freesasa.Structure and freesasa.Result objects."""
        structure = freesasa.Structure(self.topology, classifier=self.classifier)
        result = freesasa.calc(structure, parameters=self.parameters)
        return structure, result

    def atom_areas_structure(self, structure, result):
        n_atoms = structure.nAtoms()
        dtype = [("atom_name", "U10"), ("resname", "U10"), ("resid", "U10"), ("chain", "U5"), ("area", "f8")]
        arr = np.zeros(n_atoms, dtype=dtype)
        for i in range(n_atoms):
            arr[i] = (
                structure.atomName(i).strip(),
                structure.residueName(i).strip(),
                structure.residueNumber(i).strip(),
                structure.chainLabel(i).strip(),
                result.atomArea(i),
            )
        return arr

    def residue_areas_structure(self, result):
        dtype = [("chain", "U5"), ("resname", "U10"), ("resnum", "U10"), ("area", "f8")]
        rows = []
        for chain, residues in result.residueAreas().items():
            for resnum, res_area in residues.items():
                rows.append((chain, res_area.residueType.strip(), resnum, res_area.total))
        arr = np.zeros(len(rows), dtype=dtype)
        for i, row in enumerate(rows):
            arr[i] = row
        return arr

    def func_grp_structure(self, structure, result):
        func_group_sasa = {g: 0.0 for g in self.func_grp_labels}
        for i in range(structure.nAtoms()):
            resname = structure.residueName(i).strip()
            atom_name = structure.atomName(i).strip()
            atom_types = self.unique_residues.residues.get(resname)
            fg = atom_types.get(atom_name) if atom_types else None
            if fg:
                func_group_sasa[fg] += result.atomArea(i)
        dtype = [("functional_group", "U20"), ("area", "f8")]
        return np.array(list(func_group_sasa.items()), dtype=dtype)

    def total_struct_area(self, result):
        return result.totalArea()

    def run_structure(self):
        """Runs freesasa's file-based Structure/calc path once (static, no trajectory).
        Returns total area, residue areas, functional group areas, and atom areas."""
        structure, result = self.calculate_sasa_structure()

        self.structure_total_area = result.totalArea()
        self.structure_atom_area = self.atom_areas_structure(structure, result)
        self.structure_residue_area = self.residue_areas_structure(result)
        self.structure_func_grp_area = self.func_grp_structure(structure, result)

        return (
            self.structure_total_area,
            self.structure_residue_area,
            self.structure_func_grp_area,
            self.structure_atom_area,
        )

    # ---------------- calcCoord (single-frame helpers) ----------------

    def calculate_sasa_coords(self):
        """Call once per frame — uses current atomgroup.positions."""
        coords = self.atomgroup.positions.flatten().tolist()
        return freesasa.calcCoord(coords, self.radii.tolist(), self.parameters)

    def atom_areas_coords(self, result):
        n_atoms = len(self.atomgroup)
        dtype = [("atom_name", "U10"), ("resname", "U10"), ("resid", "U10"), ("chain", "U5"), ("area", "f8")]
        arr = np.zeros(n_atoms, dtype=dtype)
        for i, atom in enumerate(self.atomgroup):
            arr[i] = (atom.name, atom.resname, str(atom.resid), atom.segid, result.atomArea(i))
        return arr

    def residue_areas_coords(self, result):
        sums, names = defaultdict(float), {}
        for i, atom in enumerate(self.atomgroup):
            key = (atom.segid, atom.resid)
            sums[key] += result.atomArea(i)
            names[key] = atom.resname
        dtype = [("chain", "U5"), ("resname", "U10"), ("resnum", "U10"), ("area", "f8")]
        arr = np.zeros(len(sums), dtype=dtype)
        for i, (key, area) in enumerate(sums.items()):
            arr[i] = (key[0], names[key], str(key[1]), area)
        return arr

    def func_grp_coords(self, result):
        func_group_sasa = {g: 0.0 for g in self.func_grp_labels}
        for i, atom in enumerate(self.atomgroup):
            atom_types = self.unique_residues.residues.get(atom.resname)
            fg = atom_types.get(atom.name) if atom_types else None
            if fg:
                func_group_sasa[fg] += result.atomArea(i)
        dtype = [("functional_group", "U20"), ("area", "f8")]
        return np.array(list(func_group_sasa.items()), dtype=dtype)

    # full trajectory run (array-based, all frames)
    # one frame is a row

    def run_trajectory(self):
        n_frames = len(self.universe.trajectory)
        n_atoms = len(self.atomgroup)

        self.total_area = np.zeros(n_frames, dtype=float)
        self.residue_area = np.zeros((n_frames, len(self.residue_labels)), dtype=float)
        self.func_grp_area = np.zeros((n_frames, len(self.func_grp_labels)), dtype=float)
        self.atom_area = np.zeros((n_frames, n_atoms), dtype=float)

        for ts in self.universe.trajectory:
            result = self.calculate_sasa_coords()
            self.total_area[ts.frame] = result.totalArea()

            for i, atom in enumerate(self.atomgroup):
                area = result.atomArea(i)
                self.atom_area[ts.frame, i] = area

                res_key = (atom.segid, atom.resid, atom.resname)
                self.residue_area[ts.frame, self._residue_index[res_key]] += area

                atom_types = self.unique_residues.residues.get(atom.resname)
                fg = atom_types.get(atom.name) if atom_types else None
                if fg:
                    self.func_grp_area[ts.frame, self._func_grp_index[fg]] += area

        return self.total_area, self.residue_area, self.func_grp_area, self.atom_area

    # ---------------- sanity checks ----------------

    def print_sanity_checks(self):
        """Cheap, in-memory checks against self.atomgroup — safe to call
        regardless of trajectory (never builds a freesasa.Structure)."""

        print("MDAnalysis atomgroup atom count:", len(self.atomgroup))
        print("MDAnalysis atomgroup residue count:", len(self.atomgroup.residues))
        print("Radii array length:", len(self.radii))

        grp_counts = {}
        missing_atomgroup = []
        for atom in self.atomgroup:
            atom_types = self.unique_residues.residues.get(atom.resname)
            fg = atom_types.get(atom.name) if atom_types else None
            if fg:
                grp_counts[fg] = grp_counts.get(fg, 0) + 1
            else:
                missing_atomgroup.append((atom.resname, atom.resid, atom.name))

        total_mapped = sum(grp_counts.values())
        print("Total mapped atoms (functional groups, from atomgroup):", total_mapped)

        if total_mapped != len(self.atomgroup):
            print(
                f"  -> {len(missing_atomgroup)} unmapped atom(s) in atomgroup:"
            )
            for resname, resid, name in missing_atomgroup:
                print(f"     resname={resname:>4s} resid={resid:<6d} atom={name}")

        if self.has_trajectory:
            print("trajectory loaded — structure-based checks skipped")
            print("trajectory length/number of frames:",len(self.universe.trajectory))
            return

        structure, result = self.calculate_sasa_structure()

        print("freesasa Structure atom count:", structure.nAtoms())
        print(
            "freesasa residue count:",
            sum(len(r) for r in result.residueAreas().values()),
            )

        matched = 0
        missing_structure = []
        for i in range(structure.nAtoms()):
            resname = structure.residueName(i).strip()
            atom_name = structure.atomName(i).strip()
            fg = self.unique_residues.residues.get(resname, {}).get(atom_name)
            if fg is not None:
                matched += 1
            else:
                missing_structure.append(
                    (resname, structure.residueNumber(i).strip(), atom_name)
                )

        print("Total mapped atoms (functional groups, from structure):", matched)

        if matched != structure.nAtoms():
            print(
                f"  -> {len(missing_structure)} unmapped atom(s) in structure:"
            )
            for resname, resid, name in missing_structure:
                print(f"     resname={resname:>4s} resid={resid:<6s} atom={name}")

    # ---------------- output files ----------------

    def _header(self, file_purpose):
        lines = [
            "# " + "=" * 58,
            "# SASA \u2014 " + file_purpose,
            "# " + "-" * 58,
            "# This program computes solvent-accessible surface area (SASA)",
            "# from a molecular structure or trajectory using FreeSASA and MDAnalysis Python libraries,",
            "# broken down by atom, residue, and functional group.",
            "# " + "-" * 58,
            f"# This file: {file_purpose}",
            "# Units: Angstrom^2",
            "# Author: Dhanashri Rathod",
            f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"# Structure/Topology: {self.topology}",
        ]
        if self.has_trajectory:
            lines.append(f"# Trajectory: {self.trajectory}")
        lines += [
            f"# Classifier: {self.classifier_path or 'default'}",
            "# Algorithm: LeeRichards, probe-radius: 1.4, n-slices: 100",
            f"# Selection: {self.select}",
             f"# Residue count: {len(self.atomgroup.residues)}",
            f"# Atom count: {len(self.atomgroup)}",
            "# " + "=" * 58,
            "#",
        ]
        return "\n".join(lines) + "\n"

    def _output_paths(self):
        topo_path = Path(self.topology)
        stem = topo_path.stem

        out_dir = self.output_dir if self.output_dir is not None else topo_path.parent / "sasa"
        out_dir.mkdir(parents=True, exist_ok=True)

        return {
            "atom": out_dir / f"atom_area_{stem}.txt",
            "residue": out_dir / f"residue_area_{stem}.txt",
            "func_grp": out_dir / f"func_grp_area_{stem}.txt",
            "total": out_dir / f"total_area_{stem}.txt",
        }

    def write_output_files(self):
        """
    Writes atom-wise, residue-wise, functional-group-wise, and total
    SASA areas to text files in <topology_dir>/sasa/, named with the
    topology filename as a suffix.

    Static (pdb) run: functional-group file is a two-column list
        (functional_group, area).
    Trajectory run: functional-group file is a (n_frames x n_groups)
        table; first column is the frame number, remaining columns are
        the areas of each functional group (column order given in header).
    """
        paths = self._output_paths()

        with open(paths["atom"], "w") as f_atom, \
         open(paths["residue"], "w") as f_res, \
         open(paths["func_grp"], "w") as f_fg, \
         open(paths["total"], "w") as f_tot:

            f_atom.write(self._header("Atom-wise Areas"))
            f_atom.write("# atom_name  resname  resid  chain  area(\u00c5\u00b2)\n")

            f_res.write(self._header("Residue-wise Areas"))
            f_res.write("# chain  resname  resnum  area(\u00c5\u00b2)\n")

            f_fg.write(self._header("Functional-Group-wise Areas"))
            f_tot.write(self._header("Total Area"))

            if not self.has_trajectory:
                structure, result = self.calculate_sasa_structure()

                atom_areas = self.atom_areas_structure(structure, result)
                for row in atom_areas:
                    f_atom.write(f"{row['atom_name']:<10} {row['resname']:<8} "
                                f"{row['resid']:<8} {row['chain']:<6} "
                                f"{row['area']:}\n")

                residue_areas = self.residue_areas_structure(result)
                for row in residue_areas:
                    f_res.write(f"{row['chain']:<6} {row['resname']:<8} "
                                f"{row['resnum']:<6} {row['area']:}\n")

                f_fg.write("# functional_group  area(\u00c5\u00b2)\n")
                fg_areas = self.func_grp_structure(structure, result)
                for row in fg_areas:
                    f_fg.write(f"{row['functional_group']:<20} {row['area']:}\n")

                f_tot.write("# total_area(\u00c5\u00b2)\n")
                f_tot.write(f"{result.totalArea():}\n")

            else:
                f_tot.write("# frame  total_area(\u00c5\u00b2)\n")

                fg_names = None   # column order, taken from the first frame
                col_w = 24

                for ts in self.universe.trajectory:
                    frame = ts.frame
                    result = self.calculate_sasa_coords()

                    # ---- atom-wise (one block per frame) ----
                    f_atom.write(f"\n# Frame {frame}\n")
                    f_atom.write("# atom_name  resname  resid  chain  area(\u00c5\u00b2)\n")
                    for row in self.atom_areas_coords(result):
                        f_atom.write(f"{row['atom_name']:<10} {row['resname']:<8} "
                                     f"{row['resid']:<8} {row['chain']:<6} "
                                     f"{row['area']:}\n")

                    # ---- residue-wise (one block per frame) ----
                    f_res.write(f"\n# Frame {frame}\n")
                    f_res.write("# chain  resname  resnum  area(\u00c5\u00b2)\n")
                    for row in self.residue_areas_coords(result):
                        f_res.write(f"{row['chain']:<6} {row['resname']:<8} "
                                    f"{row['resnum']:<6} {row['area']:}\n")

                    # ---- functional-group-wise: frames x groups array ----
                    fg_areas = self.func_grp_coords(result)

                    if fg_names is None:
                        fg_names = [r["functional_group"] for r in fg_areas]
                        f_fg.write("# rows = frames, columns = functional groups; "
                                   "areas in \u00c5\u00b2\n")
                        f_fg.write(f"{'# frame':<8}" + "".join(
                            f"{n.replace(' ', '_'):>{col_w}}" for n in fg_names
                        ) + "\n")

                    f_fg.write(f"{frame:<8d}" + "".join(
                        f"{r['area']:>{col_w}}" for r in fg_areas
                    ) + "\n")

                    # ---- total ----
                    f_tot.write(f"{frame:<6} {result.totalArea():}\n")

        return paths
