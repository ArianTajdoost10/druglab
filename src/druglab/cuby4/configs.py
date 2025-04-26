from dataclasses import dataclass, field
from typing import Dict, List, Any, Tuple, Optional
from shutil import which
from collections import OrderedDict

@dataclass
class Cuby4Config:
    config: Dict[str, Any] = field(default_factory=OrderedDict)

    def get_string(self):
        return self.config_to_string(self.config, "")
    
    @classmethod
    def config_to_string(cls, config, indent):
        cstring = ""
        for k, v in config.items():
            if isinstance(v, bool):
                v = {True: "yes", False: "no"}.get(v)
            if v.__class__ not in [bool, str, dict, OrderedDict]:
                v = str(v)
            if isinstance(v, str):
                cstring += f"{indent}{k}: {v}\n"
            elif isinstance(v, dict):
                cstring += f"{indent}{k}:\n"
                indent += "  "
                cstring += cls.config_to_string(v, indent)
                indent = indent[:-2]
        return cstring

@dataclass
class Cuby4InterfaceConfig(Cuby4Config):
    interface: Optional[str] = None
    n_threads: int = 1
    
    def __post_init__(self):
        if self.interface:
            self.config["interface"] = self.interface
        self.config["cuby_threads"] = self.n_threads

@dataclass
class Cuby4CompositeInterfaceConfig(Cuby4InterfaceConfig):
    pass

@dataclass
class Cuby4JobConfig(Cuby4Config):
    job: Optional[str] = None
    
    def __post_init__(self):
        if self.job:
            self.config["job"] = self.job

@dataclass
class Cuby4MergedConfig(Cuby4Config):
    @classmethod
    def from_config_list(cls, configs):
        dicts = [c.config for c in configs]
        merged = OrderedDict()
        for d in dicts:
            for k, v in d.items():
                if k not in merged:
                    merged[k] = v
                elif isinstance(merged[k], list):
                    if isinstance(v, list): 
                        merged[k] += v
                elif isinstance(merged[k], dict):
                    if isinstance(v, dict): 
                        merged[k].update(v)
        mopac_kwds = ""
        for d in dicts:
            if "mopac_keywords" in d:
                mopac_kwds = " ".join([mopac_kwds, d["mopac_keywords"]])
        if len(mopac_kwds) > 0:
            merged["mopac_keywords"] = mopac_kwds
        return cls(merged)

@dataclass
class Cuby4EnergyJobConfig(Cuby4JobConfig):
    geometry: str = ""
    charge: int = 0
    mult: int = 1
    
    def __post_init__(self):
        super().__post_init__()
        self.job = "energy"
        self.config["job"] = "energy"
        self.config["geometry"] = self.geometry
        self.config["charge"] = self.charge
        self.config["multiplicity"] = self.mult

@dataclass
class Cuby4OptimizeJobConfig(Cuby4EnergyJobConfig):
    maxcycles: int = 50
    restart_file: Optional[str] = None
    
    def __post_init__(self):
        super().__post_init__()
        self.job = "optimize"
        self.config["job"] = "optimize"
        if self.restart_file:
            self.config["restart_file"] = self.restart_file
        self.config['maxcycles'] = self.maxcycles

@dataclass
class MOPACConfigUtils:
    keywords: List[str] = field(default_factory=list)
    
    def _update_c4mopac_keywords(self, kws):
        for kw in kws:
            if kw not in self.keywords:
                self.keywords.append(kw)
        if not hasattr(self, "config"):
            self.config = {}
        self.config["mopac_keywords"] = " ".join(self.keywords)

@dataclass
class Cuby4MOPACInterfaceConfig(Cuby4InterfaceConfig, MOPACConfigUtils):
    exe: str = "auto"
    method: str = "pm6"
    mozyme: bool = False
    setpi: bool = False
    setcharges: bool = False
    cvb: bool = False
    keywords: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        super().__post_init__()
        self.interface = "mopac"
        if self.exe == "auto":
            self.exe = which("mopac")
        #TODO delete this line
        self.exe = "/home/arian/mambaforge/envs/ASH/bin/mopac"
        if isinstance(self.mozyme, str):
            self.mozyme = {"yes": True, "no": False}.get(self.mozyme.lower())
        if isinstance(self.keywords, str):
            self.keywords = self.keywords.split()
        self._update_c4mopac_keywords(self.keywords)
        self.config["mopac_exe"] = self.exe
        self.config["method"] = self.method
        self.config["mopac_mozyme"] = self.mozyme

@dataclass
class MOPACMozymeConfig(Cuby4Config, MOPACConfigUtils):
    setpi: List[Tuple[int, int]] = field(default_factory=list)
    neg_cvb: List[Tuple[int, int]] = field(default_factory=list)
    
    def __post_init__(self):
        super().__post_init__()
        self.compile_setpi()
        self.compile_cvb()

    def compile_setpi(self):
        if len(self.setpi) == 0:
            self.config.pop("mopac_setpi", None)
            return
        setpi_str = ""
        for aid1, aid2 in self.setpi:
            setpi_str += f"\n  - {aid1};{aid2}"
        self.config["mopac_setpi"] = setpi_str
    
    def compile_cvb(self):
        if len(self.neg_cvb) == 0:
            return
        cvb_str = "CVB("
        for aid1, aid2 in self.neg_cvb:
            if cvb_str != "CVB(":
                cvb_str += ";"
            cvb_str += f"{aid1}:-{aid2}"
        cvb_str += ")"
        self._update_c4mopac_keywords([cvb_str])

@dataclass
class Cuby4AMBERInterfaceConfig(Cuby4InterfaceConfig):
    home: str = "auto"
    
    def __post_init__(self):
        super().__post_init__()
        self.interface = "amber"
        import shutil
        import os
        if self.home == "auto":
            self.home = os.path.dirname(os.path.dirname(shutil.which("sander")))
        self.config["amber_amberhome"] = self.home

@dataclass
class Cuby4QMMMInterfaceConfig(Cuby4InterfaceConfig):
    qm_config: Cuby4Config = field(default_factory=Cuby4Config)
    mm_config: Cuby4Config = field(default_factory=Cuby4Config)
    embedding: str = "mechanical"
    grad_on_point_charges: bool = False
    
    def __post_init__(self):
        super().__post_init__()
        self.interface = "qmmm"
        if isinstance(self.grad_on_point_charges, str):
            self.grad_on_point_charges = \
            {"no": False, "yes": True}.get(self.grad_on_point_charges)
        self.config["calculation_qm"] = self.qm_config.config
        self.config["calculation_mm"] = self.mm_config.config
        self.config["qmmm_embedding"] = self.embedding
        self.config["gradient_on_point_charges"] = \
            {True: "yes", False: "no"}.get(self.grad_on_point_charges)