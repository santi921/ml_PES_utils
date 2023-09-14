from tqdm import tqdm
import numpy as np
import json, argparse

from ase.io import read
from ase.constraints import FixAtoms, Hookean
from ase.ga.data import PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.constraints import FixAtoms
from ase.neighborlist import build_neighbor_list, natural_cutoffs
from ase.visualize import view


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--options", type=str, default="./bh_options.json")
    args = parser.parse_args()
    with open(args.options) as f:
        options = json.load(f)

    # model_file = options["model_file"]
    # atom_order = options["atom_order"]
    data_root = options["data_file"]
    db_file = "gadb.db"  # generalize
    # generalize
    # data_root = options["data_file"]

    atoms = read(data_root, index="0")
    cutoffs = natural_cutoffs(atoms, mult=1.5)
    neighborlist = build_neighbor_list(atoms, cutoffs=cutoffs, bothways=True)

    slab_index = [
        atom.index for atom in atoms if (atom.symbol == "Al" or atom.symbol == "O")
    ]
    cluster_index = [atom.index for atom in atoms if atom.symbol == "Pt"]

    propane_index = [atom.index for atom in atoms if atom.symbol == "C"]

    # get neighbor list using get_connectivity_

    # partition hydrogens into attached to carbon and attached to platinum
    hydrogen_index = [atom.index for atom in atoms if atom.symbol == "H"]
    hydrogen_index_carbon = []
    hydrogen_index_platinum = []
    hydrogen_index_slab = []

    for atom in hydrogen_index:
        # get neighbors of hydrogen
        indices, offsets = neighborlist.get_neighbors(
            atom,
        )
        sym_neighbors = [atoms[i].symbol for i in indices]

        if len(sym_neighbors) > 1:
            sym_neighbors.remove("H")

        if sym_neighbors[0] == "Pt":
            hydrogen_index_platinum.append(atom)
        elif sym_neighbors[0] == "Al" or sym_neighbors[0] == "O":
            hydrogen_index_slab.append(atom)
        else:
            hydrogen_index_carbon.append(atom)

    # add hydrogens_index_platinum to cluster_index
    cluster_index.extend(hydrogen_index_platinum)
    # add hydrogens_index_carbon to propane_index
    propane_index.extend(hydrogen_index_carbon)
    # add hydrogens_index_slab to slab_index
    slab_index.extend(hydrogen_index_slab)

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

    # separate the atoms in the slab and the rest
    slab = atoms[slab_index]
    cluster = atoms[cluster_index]
    propane = atoms[propane_index]

    pos = slab.get_positions()
    cell = slab.get_cell()
    p0 = np.array([0.0, 0.0, max(pos[:, 2])])
    v1 = cell[0, :] * 0.8
    v2 = cell[1, :] * 0.8
    v3 = cell[2, :] * 0.8
    v3[2] = 15.0

    # Define the composition of the atoms to optimize
    # atom_numbers = 2 * [47] + 2 * [79]
    atom_numbers = [int(atom.number) for atom in propane] + [
        int(atom.number) for atom in cluster
    ]
    # define the closest distance two atoms of a given species can be to each other
    unique_atom_types = get_all_atom_types(slab, atom_numbers)
    blmin = closest_distances_generator(
        atom_numbers=unique_atom_types, ratio_of_covalent_radii=0.7
    )

    # create the starting population
    sg = StartGenerator(
        slab=slab,
        blocks=[cluster, propane],
        blmin=blmin,
        box_to_place_in=[p0, [v1, v2, v3]],
        test_too_far=True,
        test_dist_to_slab=True,
    )

    # generate the starting population
    population_size = 50
    starting_population = []
    for i in tqdm(range(population_size)):
        starting_population.append(sg.get_new_candidate())

    d = PrepareDB(
        db_file_name=db_file, simulation_cell=slab, stoichiometry=atom_numbers
    )
    for a in starting_population:
        d.add_unrelaxed_candidate(a)


main()
