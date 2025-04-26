from dataclasses import dataclass
from typing import Dict, List, Any, Optional, Type, Tuple
from rdkit import Chem
from .general import (Cuby4GeneralEnergyCalculator,
                       generate_path_in_dir,
                         Cuby4GeneralBlock
)

from .configs import (Cuby4JobConfig, 
                      Cuby4EnergyJobConfig, 
                      Cuby4OptimizeJobConfig, 
                      Cuby4Config
)

@dataclass
class Cuby4MOPACEnergyCalculator(Cuby4GeneralEnergyCalculator):
    name: str = "c4mopacspc"
    display_name: str = "Cuby4-MOPAC Single Point Calculator"
    
    def generate_job_config(self, input_dict):
        rdkit_mol = input_dict["molecules"]
        geometry_file = generate_path_in_dir(6,
                                              self.work_dir, "_c4mecinp.pdb")
        Chem.MolToPDBFile(rdkit_mol, geometry_file)
        charge = Chem.GetFormalCharge(rdkit_mol)
        jc = Cuby4EnergyJobConfig(geometry=geometry_file,
                                   charge=charge, mult=1)
        jc.config["interface"] = "mopac"
        if hasattr(self.interface_config, "mozyme") and \
            self.interface_config.mozyme:
            pass
        return jc

@dataclass
class Cuby4MOPACEnergyOptimizer(Cuby4GeneralBlock):
    name: str = "c4mopacopt"
    display_name: str = "Cuby4-MOPAC Energy Optimizer"
    max_cycles: int = 10
    
    def generate_job_config(self, input_dict):
        rdkit_mol = input_dict["molecules"]
        geometry_file = generate_path_in_dir(6, self.work_dir, "_c4meoinp.pdb")
        restart_file = generate_path_in_dir(6, self.work_dir, "_c4meoout.pdb")
        Chem.MolToPDBFile(rdkit_mol, geometry_file)
        charge = Chem.GetFormalCharge(rdkit_mol)
        jc = Cuby4OptimizeJobConfig(
            geometry=geometry_file,
            charge=charge,
            mult=1,
            restart_file=restart_file,
            maxcycles=self.max_cycles
        )
        jc.config["interface"] = "mopac"
        return jc
    
    def generate_out_dict(self, input_dict, cuby4_out, full_config):
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])
        optimized_mol = Chem.MolFromPDBFile(full_config.config["restart_file"])
        return {
            "molecules": optimized_mol,
            "energy": energy
        }