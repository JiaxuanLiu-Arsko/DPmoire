import numpy as np
import os,sys,re
from ase.io.vasp import read_vasp, write_vasp
from ase import Atoms
from ase.build import make_supercell, sort
from .config import Config
from .structure_handler import StructureHandler
import numpy as np
import os

class EnvironmentHandler:
    struct_handler = None
    k_generator = None
    input_dir = None
    POTCAR_dir = None
    n_secs = None
    #params to be found:
    ENCUT = None
    RCUT1 = None
    RCUT2 = None
    elements = None
    top_elements = None
    bot_elements = None
    lat_vec = None
    sc = None
    def __init__(self, config:Config|dict, k_generator=None):
        if not isinstance(config, (Config, dict)):
            raise Exception("Unknown config type!")
        if k_generator is not None:
            self.k_generator = k_generator
        else:
            self.k_generator = self.gen_kpoints
        self.n_secs = config["n_sectors"]
        self.input_dir = config["input_dir"]
        self.POTCAR_dir = config["POTCAR_dir"]
        self.struct_handler = StructureHandler(config=config)
        self.lat_vec = self.struct_handler.new_struct.get_cell().array
        self.elements = self.get_elements(self.struct_handler.new_struct)
        self.top_elements = self.get_elements(self.struct_handler.top_atoms)
        self.bot_elements = self.get_elements(self.struct_handler.bot_atoms)
        self.sc = config["sc"]
        ens = []
        for element in self.elements:
            ens.append(os.popen(f"cat {config['POTCAR_dir']}/{element}/POTCAR|grep ENMAX").readlines()[0].split()[2].split(";")[0]) 
        ensf = [float(ens[k]) for k in range(len(ens))]
        self.ENCUT = np.max(ensf)
        self.RCUT1 = self.find_RCUT()
        self.RCUT2 = self.find_RCUT()

    def find_RCUT(self):
        max_a = np.max([self.struct_handler.top_atoms.get_cell().lengths()[0], 
                        self.struct_handler.bot_atoms.get_cell().lengths()[0]])
        d = self.struct_handler.d
        rcut = np.sqrt(max_a**2 + d**2)
        return rcut*1.1

    def get_elements(self, atoms:Atoms):
        element_list = atoms.get_chemical_symbols()
        elements = [element_list[0]]
        for item in element_list:
            if item == elements[-1]:
                continue
            else:
                elements.append(item)
        return elements

    def gen_POSCAR(self, out_dir:str):
        self.struct_handler.shift_all(out_dir, True, self.sc)

    def gen_init_environment(self, out_dir:str, layer:int):
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        self.gen_KPOINTS(f"{out_dir}")
        self.gen_INCAR(f"{out_dir}", f'{self.input_dir}/init_INCAR')
        if layer == 0:
            atoms = read_vasp(f'{self.input_dir}/bot_layer.poscar')
        else:
            atoms = read_vasp(f'{self.input_dir}/top_layer.poscar')
            infile = open(f"{out_dir}/INCAR", "r")
            out_str = ""
            for lines in infile:
                if len(lines.split())<=0:
                    continue
                if lines.split()[0] == "ML_ISTART":
                    out_str += "ML_ISTART = 1\n"
                else:
                    out_str += lines
            infile.close()
            with open(f"{out_dir}/INCAR", "w") as outfile:
                outfile.write(out_str)
            os.system(f"mv {out_dir}/ML_ABN {out_dir}/ML_AB")
            os.system(f"mv {out_dir}/ML_FFN {out_dir}/ML_FF")
        self.gen_POTCAR(self.get_elements(atoms), f"{out_dir}")
        atoms_sc = sort(make_supercell(prim=atoms, P=[[self.sc, 0, 0], [0, self.sc, 0], [0, 0, 1]]))
        write_vasp(f"{out_dir}/POSCAR", atoms_sc)
        os.system(f"cp {self.input_dir}/vdw_kernel.bindat {out_dir}")
    
    def gen_environment(self, INCAR, out_dir:str):
        for i in range(self.n_secs):
            for j in range(self.n_secs):
                out_dir_ij = f'{out_dir}/{i}_{j}'
                self.gen_INCAR(out_dir_ij, f'{self.input_dir}/{INCAR}')
                self.gen_KPOINTS(out_dir_ij)
                self.gen_POTCAR(self.elements, out_dir_ij)
                os.system(f"cp {self.input_dir}/vdw_kernel.bindat {out_dir_ij}")
                os.system(f"cp {self.input_dir}/ML_AB {out_dir_ij}")
                os.system(f"cp {self.input_dir}/ML_FF {out_dir_ij}")

        top_dir = f'{out_dir}/top_layer'
        top_atoms_sc = sort(make_supercell(self.struct_handler.top_atoms, [[self.sc, 0, 0], [0, self.sc, 0], [0, 0, 1]]))
        if not os.path.exists(top_dir):
            os.mkdir(top_dir)
        self.gen_INCAR(top_dir, f'{self.input_dir}/MD_monolayer_INCAR')
        self.gen_KPOINTS(top_dir)
        self.gen_POTCAR(self.get_elements(top_atoms_sc), top_dir)
        os.system(f"cp {self.input_dir}/vdw_kernel.bindat {top_dir}")
        write_vasp(f"{top_dir}/POSCAR", top_atoms_sc)

        bot_dir = f'{out_dir}/bot_layer'
        bot_atoms_sc = sort(make_supercell(self.struct_handler.bot_atoms, [[self.sc, 0, 0], [0, self.sc, 0], [0, 0, 1]]))
        if not os.path.exists(bot_dir):
            os.mkdir(bot_dir)
        self.gen_INCAR(bot_dir, f'{self.input_dir}/MD_monolayer_INCAR')
        self.gen_KPOINTS(bot_dir)
        self.gen_POTCAR(self.get_elements(bot_atoms_sc), bot_dir)
        os.system(f"cp {self.input_dir}/vdw_kernel.bindat {bot_dir}")
        write_vasp(f"{bot_dir}/POSCAR", bot_atoms_sc)
    
    def gen_POTCAR(self, elements:list, out_dir:str):
        pot_str = ""
        for element in elements:
            pot_str += f"{self.POTCAR_dir}/{element}/POTCAR "
        os.system(f"cat {pot_str} > {out_dir}/POTCAR")

    def gen_KPOINTS(self, out_dir:str):
        kpoints = self.k_generator(self.lat_vec)
        outfile = open(f"{out_dir}/KPOINTS", "w")
        outfile.write("K point mesh generated by DPmoire \n 0 \n Gamma\n")
        outfile.write(str(kpoints[0])+" "+str(kpoints[1])+" "+str(kpoints[2])+"\n")
        outfile.close()
    
    def gen_INCAR(self, out_dir:str, input_file:str):
        infile = open(input_file, "r")
        out_str = ""
        for lines in infile:
            if len(lines.split())<=0:
                continue
            if lines.split()[0] == "ENCUT":
                out_str += "ENCUT = " + str(self.ENCUT*1.5) + "\n"
            elif lines.split()[0] == "ML_RCUT1":
                out_str += "ML_RCUT1 = " + str(self.RCUT1) + "\n"
            elif lines.split()[0] == "ML_RCUT2":
                out_str += "ML_RCUT2 = " + str(self.RCUT2) + "\n" 
            elif lines.split()[0] == "LANGEVIN_GAMMA":
                out_str += "LANGEVIN_GAMMA = "
                for item in self.elements:
                    out_str += "1 "
                out_str += "\n"
            else:
                out_str += lines
        infile.close()
        with open(out_dir + "/INCAR", "w") as outfile:
            outfile.write(out_str)

    def gen_kpoints(self, lat_vec):
        kpoints = []
        for vec in lat_vec[:-1]:
            len = np.sqrt(vec[0]*vec[0] + vec[1]*vec[1])*self.sc
            kpoints.append(int(40/len))
        kpoints.append(1)
        return kpoints
