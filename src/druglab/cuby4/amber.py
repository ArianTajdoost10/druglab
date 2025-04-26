from dataclasses import dataclass
from typing import Dict, List, Any, Optional, Type, Tuple
from rdkit import Chem
import parmed
import openmm
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit import Molecule
from openmm.app import ForceField
from openff.units.openmm import to_openmm as openff_unit_to_openmm

from .general import (Cuby4GeneralEnergyCalculator,
                    generate_path_in_dir, 
                    Cuby4GeneralEnergyOptimizer
)

from .configs import (Cuby4JobConfig, 
                      Cuby4EnergyJobConfig, 
                      Cuby4OptimizeJobConfig, 
                      Cuby4AMBERInterfaceConfig, 
                      Cuby4Config
)

from .interface import Cuby4Interface, Cuby4InterfaceConfig


def rdatom2parmedatom(rdatom: Chem.Atom) -> parmed.Atom:
    pdbinfo = rdatom.GetPDBResidueInfo()
    if pdbinfo is None:
        pdbinfo = Chem.AtomPDBResidueInfo()
        pdbinfo.SetName(rdatom.GetSymbol()+str(rdatom.GetIdx()))

    return parmed.Atom(
        name=pdbinfo.GetName(),
        formal_charge=rdatom.GetFormalCharge(),
        mass=rdatom.GetMass(),
        atomic_number=rdatom.GetAtomicNum(),
    )

def rdkit2parmed(rdmol: Chem.Mol) -> parmed.Structure:
    stmol = parmed.Structure()
    coords = rdmol.GetConformer().GetPositions()

    resnum2count = dict()
    for i, rdatom in enumerate(rdmol.GetAtoms()):
        rdatom: Chem.Atom
        statom = rdatom2parmedatom(rdatom)
        pdbinfo = rdatom.GetPDBResidueInfo()
        if pdbinfo is None:
            stmol.add_atom(statom, resname="UNK", resnum=1)
        else:
            resnum = pdbinfo.GetResidueNumber()
            if resnum not in resnum2count:
                resnum2count[resnum] = len(resnum2count)+1
            resnum = resnum2count[resnum]
            stmol.add_atom(statom, 
                            pdbinfo.GetResidueName(), 
                            resnum)

    stmol.coordinates = coords
    
    for i, rdbond in enumerate(rdmol.GetBonds()):
        rdbond: Chem.Bond
        at1, at2 = rdbond.GetBeginAtomIdx(), rdbond.GetEndAtomIdx()
        stbond = parmed.Bond(stmol.atoms[at1], stmol.atoms[at2], 
                            order=rdbond.GetBondTypeAsDouble())
        stmol.bonds.append(stbond)

    res_number_list = []
    for res in stmol.residues:
        res: parmed.Residue
        res_num = res.number
        if res_num not in res_number_list:
            res_number_list.append(res.number)
        else:
            atoms = res.atoms
            [stmol.atoms.remove(atom) for atom in atoms]
            [stmol.residues[res_num-1].add_atom(atom) for atom in atoms]
            
    return stmol

def parametrizer(mol, stmol):
    molecule = Molecule.from_rdkit(mol, allow_undefined_stereo=True,
                                hydrogens_are_explicit=False)
    molecule.assign_partial_charges('mmff94')
    offtop = molecule.to_topology()
    offpos = molecule.conformers[0]
    ommtop = offtop.to_openmm()
    ommpos = openff_unit_to_openmm(offpos)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
    forcefield.registerTemplateGenerator(smirnoff.generator)

    ommsys = forcefield.createSystem(ommtop)
    parammed_struct = parmed.openmm.load_topology(stmol.topology, 
                                                    ommsys, 
                                                    stmol.coordinates)
    return parammed_struct


@dataclass
class Cuby4AMBEREnergyCalculator(Cuby4GeneralEnergyCalculator):
    name: str = "c4amberspc"
    display_name: str = "Cuby4-AMBER Single Point Calculator"
    
    def generate_job_config(self, input_dict: Dict[str, Any])\
        -> Cuby4JobConfig:
        mol = input_dict["molecules"]
        geometry_file = generate_path_in_dir(6, self.work_dir, ".pdb")
        prmtop_file = generate_path_in_dir(6, self.work_dir, ".parm7")
        
        Chem.MolToPDBFile(mol, geometry_file)
        
        stmol = rdkit2parmed(mol)
        struct = parametrizer(mol, stmol)
        struct.save(prmtop_file)
        
        config = Cuby4EnergyJobConfig(geometry=geometry_file)
        config.config["amber_top_file"] = prmtop_file
        config.config["interface"] = "amber"
        config.config["job"] = "energy"
        return config

@dataclass
class Cuby4AMBEREnergyOptimizer(Cuby4GeneralEnergyOptimizer):
    name: str = "c4amberopt"
    display_name: str = "Cuby4-AMBER Energy Optimizer"
    
    def generate_job_config(self, input_dict: Dict[str, Any])\
          -> Cuby4JobConfig:
        mol = input_dict["molecules"]
        geometry_file = generate_path_in_dir(6, self.work_dir, ".pdb")
        prmtop_file = generate_path_in_dir(6, self.work_dir, ".parm7")
        restart_file = generate_path_in_dir(6, self.work_dir, ".pdb")

        Chem.MolToPDBFile(mol, geometry_file)

        stmol = rdkit2parmed(mol)
        struct = parametrizer(mol, stmol)
        struct.save(prmtop_file)
        
        config = Cuby4OptimizeJobConfig(
            geometry=geometry_file,
            restart_file=restart_file
        )
        config.config["amber_top_file"] = prmtop_file
        config.config["interface"] = "amber"
        config.config["job"] = "optimize"
        return config
    
@dataclass
class Cuby4AMBERMultiComponentEnergyOptimizer(Cuby4GeneralEnergyOptimizer):
    name: str = "c4ambermcopt"
    display_name: str = "Cuby4-AMBER Multi Component Energy Optimizer"
    num_components: int = 1
    
    def generate_job_config(self, input_dict: Dict[str, Any])\
          -> Cuby4JobConfig:
        all_mols = [input_dict[f"molecule_{i+1}"] 
                    for i in range(self.num_components)]
        
        combined_mol = Chem.Mol(all_mols[0])
        for mol in all_mols[1:]:
            combined_mol = Chem.CombineMols(combined_mol, mol)
            
        geometry_file = generate_path_in_dir(6, self.work_dir, ".pdb")
        prmtop_file = generate_path_in_dir(6, self.work_dir, ".parm7")
        restart_file = generate_path_in_dir(6, self.work_dir, ".pdb")
        
        Chem.MolToPDBFile(combined_mol, geometry_file)
        
        stmol = rdkit2parmed(mol)
        struct = parametrizer(mol, stmol)
        struct.save(prmtop_file)
        
        config = Cuby4OptimizeJobConfig(
            geometry=geometry_file,
            restart_file=restart_file
        )
        config.config["amber_top_file"] = prmtop_file
        config.config["interface"] = "amber"
        config.config["job"] = "optimize"
        return config
    
    def generate_out_dict(self, 
                          input_dict: Dict[str, Any], 
                          cuby4_out: str, 
                          full_config: Cuby4Config) -> Dict[str, Any]:
        energy = float(cuby4_out.split("Energy:")[-1].split()[0])
        optimized_mol = Chem.MolFromPDBFile(full_config.config["restart_file"])
        
        output_dict = {}
        for i in range(self.num_components):
            output_dict[f"molecule_{i+1}"] = input_dict[f"molecule_{i+1}"]
        return output_dict