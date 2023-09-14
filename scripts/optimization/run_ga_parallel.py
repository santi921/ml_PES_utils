from random import random
import json
import argparse
from ase.io import write
import time
from ase.ga.data import DataConnection
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.parallellocalrun import ParallelLocalRun
from nequip.ase import NequIPCalculator


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--options", type=str, default="./bh_options.json")
    args = parser.parse_args()
    with open(args.options) as f:
        options = json.load(f)

    population_size = options["population_size"]
    mutation_probability = options["mutation_probability"]
    n_to_test = options["n_to_test"]
    model_file = options["model_file"]
    atom_order = options["atom_order"]

    model_str = model_file
    config = {
        "chemical_symbols": atom_order,
    }
    nequip_calc = NequIPCalculator.from_deployed_model(
        model_str,
        set_global_options=True,
        species_to_type_name={s: s for s in config["chemical_symbols"]},
    )

    # Initialize the different components of the GA
    da = DataConnection("gadb.db")
    tmp_folder = "tmp_folder/"

    # An extra object is needed to handle the parallel execution
    parallel_local_run = ParallelLocalRun(
        data_connection=da,
        tmp_folder=tmp_folder,
        n_simul=4,
        calc_script="./run_ga_single.py",
    )

    atom_numbers_to_optimize = da.get_atom_numbers_to_optimize()
    n_to_optimize = len(atom_numbers_to_optimize)
    slab = da.get_slab()
    all_atom_types = get_all_atom_types(slab, atom_numbers_to_optimize)
    blmin = closest_distances_generator(all_atom_types, ratio_of_covalent_radii=0.7)

    comp = InteratomicDistanceComparator(
        n_top=n_to_optimize,
        pair_cor_cum_diff=0.015,
        pair_cor_max=0.7,
        dE=0.02,
        mic=False,
    )
    pairing = CutAndSplicePairing(slab, n_to_optimize, blmin)
    mutations = OperationSelector(
        [1.0, 1.0, 1.0],
        [
            MirrorMutation(blmin, n_to_optimize),
            RattleMutation(blmin, n_to_optimize),
            PermutationMutation(n_to_optimize),
        ],
    )

    # Relax all unrelaxed structures (e.g. the starting population)
    while da.get_number_of_unrelaxed_candidates() > 0:
        a = da.get_an_unrelaxed_candidate()
        parallel_local_run.relax(a)

    # Wait until the starting population is relaxed
    while parallel_local_run.get_number_of_jobs_running() > 0:
        time.sleep(5.0)

    # create the population
    population = Population(
        data_connection=da, population_size=population_size, comparator=comp
    )

    # test n_to_test new candidates
    for i in range(n_to_test):
        print("Now starting configuration number {0}".format(i))
        a1, a2 = population.get_two_candidates()
        a3, desc = pairing.get_new_individual([a1, a2])
        if a3 is None:
            continue
        da.add_unrelaxed_candidate(a3, description=desc)

        # Check if we want to do a mutation
        if random() < mutation_probability:
            a3_mut, desc = mutations.get_new_individual([a3])
            if a3_mut is not None:
                da.add_unrelaxed_step(a3_mut, desc)
                a3 = a3_mut

        # Relax the new candidate
        parallel_local_run.relax(a3)
        population.update()

    # Wait until the last candidates are relaxed
    while parallel_local_run.get_number_of_jobs_running() > 0:
        time.sleep(5.0)

    write("all_candidates.traj", da.get_all_relaxed_candidates())


main()
