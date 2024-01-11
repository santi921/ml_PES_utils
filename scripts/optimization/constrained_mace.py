import json, argparse
import torch
from ase.io import read
from ase.constraints import FixAtoms, Hookean
from ase.optimize.minimahopping import MinimaHopping
from mace.calculators import mace_mp
#from nequip.ase import NequIPCalculator


def main():
    # data_root = "/home/santiagovargas/dev/ml_PES_utils/data/bh_09132023/input.traj"

    parser = argparse.ArgumentParser()
    parser.add_argument("--options", type=str, default="./bh_options.json")
    args = parser.parse_args()
    with open(args.options) as f:
        options = json.load(f)

    #model_file = options["model_file"]
    atom_order = options["atom_order"]
    data_root = options["data_file"]

    atoms = read(data_root, index="0")

    slab_index = [
        atom.index for atom in atoms if (atom.symbol == "Al" or atom.symbol == "O")
    ]
    non_slab_index = [
        atom.index for atom in atoms if (atom.symbol != "Al" and atom.symbol != "O")
    ]

    constraints = [
        FixAtoms(
            indices=[
                atom.index
                for atom in atoms
                if (atom.symbol == "Al" or atom.symbol == "O")
            ]
        )
    ]

    dict_bonds = {
        "C-Pt": 1.9,
        "Pt-Pt": 2.869,
        "C-C": 1.54,
        "C-O": 1.43,
        "Pt-O": 2.0,
        "O-O": 1.2,
        "Al-Al": 2.39,
        "H-H": 0.74,
        "H-Pt": 1.89,
        "H-C": 1.09,
        "Al-O": 1.87,
        "Al-Pt": 2.39,
        "Al-H": 1.66,
        "Al-C": 2.13,
        "O-H": 0.96,
    }

    dict_k = {
        "C-Pt": 3,
        "Pt-Pt": 2,
        "C-C": 6,
        "C-O": 5,
        "Pt-O": 2,
        "O-O": 5,
        "Al-Al": 3,
        "H-H": 5,
        "H-Pt": 4,
        "H-C": 7,
        "Al-O": 4,
        "Al-Pt": 2,
        "Al-H": 3,
        "Al-C": 3,
        "O-H": 5,
    }

    for i in non_slab_index:
        for j in non_slab_index:
            if i != j:
                # get the elemnts of the atoms
                elem1 = atoms[i].symbol
                elem2 = atoms[j].symbol
                key_dict = elem1 + "-" + elem2
                key_rev = elem2 + "-" + elem1
                if key_rev in dict_bonds.keys():
                    key_dict = key_rev

                constraints.append(
                    Hookean(a1=i, a2=j, rt=dict_bonds[key_dict], k=dict_k[key_dict])
                )

    
    config = {
        "chemical_symbols": atom_order,
    }
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    #nequip_calc = NequIPCalculator.from_deployed_model(
    #    model_str,
    #    set_global_options=True,
    #    device=device,
    #    species_to_type_name={s: s for s in config["chemical_symbols"]},
    #)
    #atoms.set_calculator(nequip_calc)
    calc = mace_mp(model="large", device="cuda")
    atoms.set_calculator(calc)
    hop = MinimaHopping(atoms, Ediff0=1.0, T0=4000.0)
    hop(totalsteps=100)


main()
