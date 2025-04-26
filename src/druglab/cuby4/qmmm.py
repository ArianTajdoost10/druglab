from dataclasses import dataclass
from typing import Dict, List, Any, Optional, Type, Tuple
from rdkit import Chem
from .amber import rdkit2parmed, rdatom2parmedatom, parametrizer
import parmed
from .general import (Cuby4GeneralBlock, 
                      Cuby4GeneralEnergyCalculator, 
                      generate_path_in_dir, 
                      Cuby4GeneralEnergyOptimizer
)

from .configs import (Cuby4Config, 
                      Cuby4JobConfig, 
                      Cuby4OptimizeJobConfig, 
                      Cuby4EnergyJobConfig
)

@dataclass
class Cuby4QMMMEnergyCalculator(Cuby4GeneralEnergyCalculator):
    name: str = "c4qmmmspc"
    display_name: str = "Cuby4-QMMM Single Point Calculator"
    
    def generate_job_config(self, input_dict):
        qm_mol = input_dict["qm_region"]
        mm_mol = input_dict["nonqm_region"]
        geometry_file = generate_path_in_dir(6, self.work_dir, ".pdb")
        combined_mol = Chem.CombineMols(qm_mol, mm_mol)
        Chem.MolToPDBFile(combined_mol, geometry_file)
        parm7_whole = generate_path_in_dir(6, self.work_dir, ".parm7")
        parm7_qm = generate_path_in_dir(6, self.work_dir, ".parm7")

        stmol = rdkit2parmed(combined_mol)
        struct: parmed.Structure = parametrizer(combined_mol, stmol)
        struct.save(parm7_whole)

        stmol = rdkit2parmed(qm_mol)
        struct: parmed.Structure = parametrizer(qm_mol, stmol)
        struct.save(parm7_qm)

        qm_size = qm_mol.GetNumAtoms()
        charge_whole = Chem.GetFormalCharge(combined_mol)
        charge_qm = Chem.GetFormalCharge(qm_mol)
        config = Cuby4EnergyJobConfig(geometry=geometry_file,
                                       charge=charge_whole)
        config.config["qmmm_core"] = f"1-{qm_size}"
        config.config["calculation_qm"] = {}
        config.config["calculation_qm"]["interface"] = "mopac"
        config.config["calculation_qm"]["charge"] = charge_qm
        config.config["calculation_mm"] = {}
        config.config["calculation_mm"]["interface"] = "amber"
        config.config["calculation_mm"]["amber_top_file"] = parm7_whole
        config.config["calculation_qmregion_mm"] = {"amber_top_file": parm7_qm}
        config.config["interface"] = "qmmm"
        return config

@dataclass
class Cuby4QMMMEnergyOptimizer(Cuby4GeneralEnergyOptimizer):
    name: str = "c4qmmmspcop"
    display_name: str = "Cuby4-QMMM Energy Optimizer"
    maxcycles: int = 10
    
    def generate_job_config(self, input_dict):
        qm_mol = input_dict["qm_region"]
        mm_mol = input_dict["nonqm_region"]
        geometry_file = generate_path_in_dir(6, self.work_dir, ".pdb")
        restart_file = generate_path_in_dir(6, self.work_dir, ".pdb")
        combined_mol = Chem.CombineMols(qm_mol, mm_mol)
        Chem.MolToPDBFile(combined_mol, geometry_file)
        parm7_whole = generate_path_in_dir(6, self.work_dir, ".parm7")
        parm7_qm = generate_path_in_dir(6, self.work_dir, ".parm7")
        with open(parm7_whole, "w") as f:
            f.write("# Placeholder topology file for whole system\n")
        with open(parm7_qm, "w") as f:
            f.write("# Placeholder topology file for QM region\n")
        qm_size = qm_mol.GetNumAtoms()
        charge_whole = Chem.GetFormalCharge(combined_mol)
        charge_qm = Chem.GetFormalCharge(qm_mol)
        config = Cuby4OptimizeJobConfig(
            geometry=geometry_file,
            charge=charge_whole,
            maxcycles=self.maxcycles,
            restart_file=restart_file
        )
        config.config["qmmm_core"] = f"1-{qm_size}"
        config.config["calculation_qm"] = {}
        config.config["calculation_qm"]["charge"] = charge_qm
        config.config["calculation_qm"]["interface"] = "mopac"
        config.config["calculation_mm"] = {}
        config.config["calculation_mm"]["interface"] = "amber"
        config.config["calculation_mm"]["amber_top_file"] = parm7_whole
        config.config["calculation_qmregion_mm"] = {"amber_top_file": parm7_qm}
        config.config["interface"] = "qmmm"
        return config
    
    def generate_out_dict(self, input_dict, cuby4_out, full_config):
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])
        qm_size = input_dict["qm_region"].GetNumAtoms()
        qm_region = input_dict["qm_region"]
        nonqm_region = input_dict["nonqm_region"]
        return {
            "energy": energy,
            "qm_region": qm_region,
            "nonqm_region": nonqm_region
        }