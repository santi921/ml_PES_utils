# Convert coded species to 0, 1, 2, etc.
coded_species = [species_map[spec] for spec in atoms.numbers]


# read in atoms via ase



import numpy as np
import matplotlib.pyplot as plt


def eval(sgp_calc, test_structures):
    energy_correlation = []
    force_correlation = []
    # sgp_calc = SGP_Calculator(sparse_gp)
    for atoms in test_structures:
        sgp_calc.calculate(atoms)
        atoms.calc.calculate(atoms)
        energy_correlation.append(
            (atoms.calc.results["energy"], sgp_calc.results["energy"])
        )
        for force_calc, force_sgp in zip(
            atoms.calc.results["forces"], sgp_calc.results["forces"]
        ):
            force_correlation.append((force_calc, force_sgp))
    energy_correlation = np.array(energy_correlation)
    force_correlation = np.array(force_correlation)

    fig, axs = plt.subplots(1, 2, figsize=(7, 7))
    for ax, correlation, title in zip(
        axs, [energy_correlation, force_correlation], ["energy", "force"]
    ):
        ax.scatter(correlation[:, 0], correlation[:, 1])
        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()]),
            np.max([ax.get_xlim(), ax.get_ylim()]),
        ]

        ax.plot(lims, lims, "k-", alpha=0.75, zorder=0)
        ax.set_aspect("equal")
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_title(title)
