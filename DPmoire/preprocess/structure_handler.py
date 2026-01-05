import numpy as np
import os
from .config import Config
from ase.io.vasp import read_vasp, write_vasp
from ase.constraints import FixedLine
from ase.build import make_supercell, sort, stack
from ase import Atoms
import copy
from ._find_homo_twist import search_twist, adjust_atoms_d
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.ase import AseAtomsAdaptor


class StructureHandler:
    """
    Handle crystal structures during the process.
    """

    in_dir = None
    out_dir = None  # Directory to store shifted structures.
    work_dir = None
    n_secs_layer2 = None
    n_secs_layer3 = None
    layer_atoms = None
    layer_indices = None
    layer_centers_frac = None
    new_struct = None
    d = None

    def __init__(self, config: Config | dict = None):
        if isinstance(config, (Config, dict)):
            self.n_secs_layer2 = config["n_sectors_layer2"]
            self.n_secs_layer3 = config["n_sectors_layer3"]
            self.in_dir = config["input_dir"]
            self.work_dir = config["work_dir"]
            self.d = config["d"]
        else:
            raise Exception(f"Unknown type of Conifg:{type(config)}")
        self.read_all_layers(self.in_dir)
        self.new_struct, self.layer_indices = self.build_new_struct(d=self.d)

    def read_atoms(self, in_file: str):
        atoms = read_vasp(in_file)
        return atoms

    def read_all_layers(self, in_dir: str):
        self.layer_atoms = [
            self.read_atoms(f"{in_dir}/layer1.poscar"),
            self.read_atoms(f"{in_dir}/layer2.poscar"),
            self.read_atoms(f"{in_dir}/layer3.poscar"),
        ]
        # Keep backwards-compatible aliases for other modules.
        self.bot_atoms = self.layer_atoms[0]
        self.mid_atoms = self.layer_atoms[1]
        self.top_atoms = self.layer_atoms[2]

    def _fractional_positions(self, atoms: Atoms):
        frac_mat = np.linalg.inv(atoms.get_cell().array)
        return np.dot(atoms.get_positions(), frac_mat)

    def _assign_layer_indices(self, atoms: Atoms, centers_frac: list[float]):
        frac_pos = self._fractional_positions(atoms)
        indices = [[] for _ in centers_frac]
        for idx, pos in enumerate(frac_pos):
            z = pos[2] - np.floor(pos[2])
            dists = [
                min(abs(z - center), 1 - abs(z - center)) for center in centers_frac
            ]
            layer_idx = int(np.argmin(dists))
            indices[layer_idx].append(idx)
        return indices

    def _generate_all_stackings(self):
        stackings = []
        seen = set()
        for i2 in range(self.n_secs_layer2):
            for j2 in range(self.n_secs_layer2):
                vec = (i2, j2, 0, 0)
                if vec not in seen:
                    stackings.append([i2, j2, 0, 0])
                    seen.add(vec)
        for i3 in range(self.n_secs_layer3):
            for j3 in range(self.n_secs_layer3):
                vec = (0, 0, i3, j3)
                if vec not in seen:
                    stackings.append([0, 0, i3, j3])
                    seen.add(vec)
        return stackings

    def find_sym_reduced_stackings(self, prec: float = 0.0001):
        _ = prec  # kept for API compatibility
        adaptor = AseAtomsAdaptor()
        matcher = StructureMatcher(ltol = prec, stol=prec, angle_tol=prec)
        unique_structs = []
        unique_stackings = []
        for stacking in self._generate_all_stackings():
            atoms = self._shift_primitive(*stacking)
            structure = adaptor.get_structure(atoms)
            matched = False
            for ref in unique_structs:
                if matcher.fit(ref, structure):
                    matched = True
                    break
            if not matched:
                unique_structs.append(structure)
                unique_stackings.append(stacking)
        stcks = np.array(unique_stackings)
        np.savetxt(f"{self.work_dir}/sym_reduced_stackings.txt", stcks, fmt="%d")
        return stcks

    def build_new_struct(self, d: float):
        """
        Build a trilayer structure combining all three layers.
        """
        n_layers = len(self.layer_atoms)
        ref_cell = self.layer_atoms[0].get_cell().array
        layer_lengths = [atoms.get_cell().lengths() for atoms in self.layer_atoms]
        new_cell_mat = np.zeros_like(ref_cell)
        for axis in range(3):
            ref_vec = ref_cell[axis]
            ref_len = np.linalg.norm(ref_vec)
            if ref_len < 1e-8:
                raise ValueError("Invalid lattice vector length encountered.")
            target_len = np.mean([lengths[axis] for lengths in layer_lengths])
            new_cell_mat[axis] = ref_vec / ref_len * target_len
        c_len = np.linalg.norm(new_cell_mat[2])
        spacing_frac = d / c_len
        centers_frac = [
            0.5 + (idx - (n_layers - 1) / 2) * spacing_frac for idx in range(n_layers)
        ]
        new_pos = []
        new_symbols = []
        for layer_idx, atoms in enumerate(self.layer_atoms):
            frac = self._fractional_positions(atoms)
            frac[:, 2] += centers_frac[layer_idx] - np.mean(frac[:, 2])
            cart = np.dot(frac, new_cell_mat)
            new_pos.extend(cart.tolist())
            new_symbols.extend(atoms.get_chemical_symbols())
        atoms = sort(
            Atoms(
                positions=new_pos,
                symbols=new_symbols,
                cell=new_cell_mat,
                pbc=[True, True, True],
            )
        )
        self.layer_centers_frac = centers_frac
        layer_indices = self._assign_layer_indices(atoms, centers_frac)
        return atoms, layer_indices

    def _shift_primitive(self, i2: int, j2: int, i3: int, j3: int):
        atoms = copy.deepcopy(self.new_struct)
        pos = atoms.get_positions()
        cell = atoms.get_cell().array
        delta2 = (
            i2 / self.n_secs_layer2 * cell[0] + j2 / self.n_secs_layer2 * cell[1]
        )
        delta3 = (
            i3 / self.n_secs_layer3 * cell[0] + j3 / self.n_secs_layer3 * cell[1]
        )
        for idx in self.layer_indices[1]:
            pos[idx] += delta2
        for idx in self.layer_indices[2]:
            pos[idx] += delta3
        atoms.set_positions(pos)
        return atoms

    def shift_atoms(
        self, i2: int, j2: int, i3: int, j3: int, c_constrain: bool = True, sc: int = 2
    ):
        """
        Return atoms where the 2nd layer is shifted by (i2/n2 * a1 + j2/n2 * a2) and
        the 3rd layer is shifted by (i3/n3 * a1 + j3/n3 * a2).
        """
        atoms = self._shift_primitive(i2, j2, i3, j3)
        atoms_sc = sort(
            make_supercell(prim=atoms, P=[[sc, 0, 0], [0, sc, 0], [0, 0, 1]])
        )
        layer_indices_sc = self._assign_layer_indices(
            atoms_sc, self.layer_centers_frac
        )
        if c_constrain and len(layer_indices_sc[0]) > 0 and len(layer_indices_sc[1]) > 0 and len(layer_indices_sc[2]) > 0:
            cons = FixedLine(
                [layer_indices_sc[0][0], layer_indices_sc[1][0], layer_indices_sc[2][0]],
                direction=atoms_sc.cell.array[2] / atoms_sc.cell.lengths()[2],
            )
            atoms_sc.set_constraint(cons)
        return atoms_sc

    def shift(self, stacking, out_dir: str, c_constrain: bool = True, sc: int = 2):
        """
        Write POSCAR at out_dir/i2_j2_i3_j3/ of shifted structures.
        """
        i2, j2, i3, j3 = map(int, stacking)
        shift_dir = f"{out_dir}/{i2}_{j2}_{i3}_{j3}/"
        if not os.path.exists(shift_dir):
            os.makedirs(shift_dir)
        atoms_sc = self.shift_atoms(i2, j2, i3, j3, c_constrain, sc)
        write_vasp(f"{shift_dir}/POSCAR", atoms=atoms_sc)

    def shift_all(
        self,
        out_dir: str,
        c_constrain: bool = True,
        sc: int = 2,
        stackings=None,
    ):
        if stackings is None:
            stackings = self._generate_all_stackings()
        for stck in stackings:
            self.shift(stck, out_dir, c_constrain=c_constrain, sc=sc)

    def make_twist_struct(self, N_min, N_max, out_dir: str):
        top_atoms = copy.deepcopy(self.top_atoms)
        bot_atoms = copy.deepcopy(self.bot_atoms)
        top_atoms, bot_atoms = adjust_atoms_d(top_atoms, bot_atoms, self.d)
        angle_list, mat_list = search_twist(N_min, N_max)
        out_atoms_list = []
        for idx, mat in enumerate(mat_list):
            top_sc = make_supercell(top_atoms, P=mat[0])
            bot_sc = make_supercell(bot_atoms, P=mat[1])
            out_atoms = stack(bot_sc, top_sc, maxstrain=None, reorder=True)
            out_atoms_list.append(out_atoms)
            if not os.path.exists(f"{out_dir}/{angle_list[idx]}/"):
                os.mkdir(f"{out_dir}/{angle_list[idx]}/")
            write_vasp(f"{out_dir}/{angle_list[idx]}/POSCAR", out_atoms)
        return angle_list, out_atoms_list
