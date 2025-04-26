from dataclasses import dataclass
import pathlib, shutil, subprocess, os, glob
import string
import random
from typing import List, Tuple, Dict, Any, Optional, Union
from rdkit import Chem
from .configs import Cuby4InterfaceConfig, Cuby4Config, Cuby4MergedConfig

def generate_name_in_dir(length, dir_path, extension):
    chars = string.ascii_lowercase + string.digits
    while True:
        name = ''.join(random.choice(chars) for _ in range(length))
        full_path = os.path.join(dir_path, f"{name}{extension}")
        if not os.path.exists(full_path):
            return full_path

@dataclass
class Cuby4Interface:
    interface_config: Cuby4InterfaceConfig
    cuby4_exe: str = "auto"
    work_dir: str = "."
    
    def __post_init__(self):
        self.work_dir = str(pathlib.Path(self.work_dir).absolute())
        self.path = shutil.which("cuby4") if \
            self.cuby4_exe == "auto" else \
                str(pathlib.Path(self.cuby4_exe).absolute())
    
    def run(self, config, job_name=None):
        if job_name is None:
            job_name = generate_name_in_dir(6, self.work_dir, ".yml")
        config_path = os.path.join(self.work_dir, f"{job_name}_conf.yml")
        config_path = str(pathlib.Path(config_path).absolute())
        config_final = Cuby4MergedConfig.from_config_list([self.interface_config,
                                                            config])
        config_str = config_final.get_string()
        with open(config_path, "w") as f:
            f.write(config_str)
        prev_path = pathlib.Path(os.curdir).absolute()
        os.chdir(self.work_dir)
        for filename in glob.glob("job_*"):
            shutil.rmtree(filename) 
        output = subprocess.run([self.path, config_path], capture_output=True)
        os.chdir(prev_path)
        if output.returncode != 0:
            print(output)
            raise ValueError("Cuby4 execution failed")
        return output.stdout.decode()