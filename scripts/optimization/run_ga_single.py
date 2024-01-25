from ase.optimize import BFGS
from ase.io import read, write
from ase.calculators.emt import EMT
from ase.ga.relax_attaches import VariansBreak
import sys
import json
import torch

from nequip.ase import NequIPCalculator


def main():
    fname = sys.argv[1]

    with open("./options.json") as f:
        options = json.load(f)

    model_file = options["model_file"]
    atom_order = options["atom_order"]

    model_str = model_file
    config = {
        "chemical_symbols": atom_order,
    }
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    nequip_calc = NequIPCalculator.from_deployed_model(
        model_str,
        set_global_options=True,
        device=device,
        species_to_type_name={s: s for s in config["chemical_symbols"]},
    )

    print("Now relaxing {0}".format(fname))
    a = read(fname)
    a.set_calculator(nequip_calc)
    dyn = BFGS(a, trajectory=None, logfile=None)
    vb = VariansBreak(a, dyn)
    dyn.attach(vb.write)
    dyn.run(fmax=0.05)

    a.info["key_value_pairs"]["raw_score"] = -a.get_potential_energy()

    write(fname[:-5] + "_done.traj", a)

    print("Done relaxing {0}".format(fname))


main()
