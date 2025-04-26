from dataclasses import dataclass, field
import os
import pathlib
from typing import Dict, List, Any, Optional, Type, Tuple, Callable
from rdkit import Chem
import string, random
from .configs import Cuby4Config, Cuby4MergedConfig, Cuby4JobConfig
from .interface import Cuby4Interface, Cuby4InterfaceConfig

def generate_path_in_dir(length, dir_path, extension):
    chars = string.ascii_lowercase + string.digits
    while True:
        name = ''.join(random.choice(chars) for _ in range(length))
        full_path = os.path.join(dir_path, f"{name}{extension}")
        if not os.path.exists(full_path):
            return full_path
        
@dataclass
class Cuby4GeneralBlock:
    interface_config: Cuby4InterfaceConfig
    work_dir: str = "."
    cuby4_exe: str = "auto"
    debug: bool = False
    save_output: bool = False
    
    def __post_init__(self):
        self.interface = Cuby4Interface(
            interface_config=self.interface_config,
            work_dir=self.work_dir,
            cuby4_exe=self.cuby4_exe
        )

    def operate(self, input_dict: Dict[str, Any]) -> Dict[str, Any]:
        jc = self.generate_job_config(input_dict)
        full_config = Cuby4MergedConfig.from_config_list([self.interface_config, 
                                                          jc])
        output = self.interface.run(full_config)
        return self.generate_out_dict(input_dict, output, full_config)

    def generate_job_config(self, input_dict: Dict[str, Any]) -> Cuby4JobConfig:
        raise NotImplementedError()
    
    def generate_out_dict(self, 
                          input_dict: Dict[str, Any], 
                          cuby4_out: str,
                          final_config: Cuby4Config) -> Dict[str, Any]:
        raise NotImplementedError()

@dataclass
class Cuby4GeneralEnergyCalculator(Cuby4GeneralBlock):
    name: str = "c4general"
    display_name: str = "Cuby4 General Energy Calculator"
    
    def generate_out_dict(self, 
                          input_dict: Dict[str, Any], 
                          cuby4_out: str, 
                          final_config: Cuby4Config) -> Dict[str, Any]:
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])
        return {"energy": energy}

@dataclass
class Cuby4GeneralEnergyOptimizer(Cuby4GeneralBlock):
    def generate_out_dict(self, 
                          input_dict: Dict[str, Any], 
                          cuby4_out: str, 
                          full_config: Cuby4Config) -> Dict[str, Any]:
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])
        optimized_mol = Chem.MolFromPDBFile(full_config.config["restart_file"])
        return {
            "molecules": optimized_mol,
            "energy": energy
        }