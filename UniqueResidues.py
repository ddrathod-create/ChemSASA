import numpy as np
import MDAnalysis as mda

from res_func_grp import BACKBONE, SIDECHAINS

class UniqueResidues:
    def __init__(self):
        self.residues = {
            resname: {**BACKBONE, **sidechain}
            for resname, sidechain in SIDECHAINS.items()
        }

    def count_func_grp_from_universe(self, universe):
        grp_counts = {}

        for atom in universe.atoms:
            atom_types = self.residues.get(atom.resname)
            if atom_types is None:
                continue
            functional_group = atom_types.get(atom.name)
            if functional_group is None:
                continue
            grp_counts[functional_group] = grp_counts.get(functional_group, 0) + 1

        dtype = [("functional_group", "U20"), ("count", "i4")]
        counts_array = np.array(list(grp_counts.items()), dtype=dtype)
        total = np.array(counts_array["count"].sum())

        return {"counts": counts_array, "total": total}

    def count_func_grp_from_structure(self, structure):
        grp_counts = {}

        for i in range(structure.nAtoms()):
            resname = structure.residueName(i).strip()
            atom_name = structure.atomName(i).strip()

            atom_types = self.residues.get(resname)
            if atom_types is None:
                continue
            functional_group = atom_types.get(atom_name)
            if functional_group is None:
                continue
            grp_counts[functional_group] = grp_counts.get(functional_group, 0) + 1

        dtype = [("functional_group", "U20"), ("count", "i4")]
        counts_array = np.array(list(grp_counts.items()), dtype=dtype)
        total = np.array(counts_array["count"].sum())

        return {"counts": counts_array, "total": total}
