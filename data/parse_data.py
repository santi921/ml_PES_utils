import os
import re
import numpy as np

# import matplotlib.pyplot as plt
# import seaborn as sns
import pandas as pd


def read_xyz(file="./gas_phase/C3H6/c3h6.xyz"):
    """Read the xyz file and return the coordinates of the atoms"""
    with open(file, "r") as f:
        lines = f.readlines()
    return lines


def parse_lines(lines):
    # find lines that start with a number
    # and split them into a list of strings
    dict_energies = {"energy": [], "force": [], "stress": []}

    ind = -1
    start_control = False

    for line in lines:
        line_split = line.split()
        # print(line_split)

        if line_split[0].isdigit():
            pass

        elif line[0:7] == "Lattice":
            if start_control:
                # print(np.mean(forces))
                dict_energies["force"].append(np.mean(forces))
            forces = []

            start_control = True

            # find the index of energy string in line
            # and take line from there to the end

            # check if energy in line
            if "energy" in line:
                ind_energy_start = line.index("energy")
                white_spaces_ind = [m.start() for m in re.finditer(" ", line)]
                ind_energy_end = white_spaces_ind[0]

                for i in white_spaces_ind:
                    if i > ind_energy_start:
                        ind_energy_end = i
                        break

                energy = float(line[ind_energy_start:ind_energy_end].split("=")[1])

                dict_energies["energy"].append(energy)

            if "stress" in line:
                ind_stress_start = line.index("stress")
                # find next two quotes after ind_stress_start
                ind_stress_str_start = line[ind_stress_start:].index('"')
                ind_stress_str_end = line[
                    ind_stress_start + ind_stress_str_start + 1 :
                ].index('"')
                stress = [
                    float(i)
                    for i in line[
                        ind_stress_str_start
                        + ind_stress_start
                        + 1 : ind_stress_start
                        + ind_stress_str_end
                    ].split(" ")
                ]

                dict_energies["stress"].append(
                    np.linalg.norm(stress)
                )  # np.array([np.linalg.norm(np.array(i)) for i in stress])

            # if "magmom" in line:
            #    ind_magmom_start = line.index("magmom")
            #    white_spaces_ind = [m.start() for m in re.finditer(' ', line)]
            #    ind_magmom_end = white_spaces_ind[0]#

            #    for i in white_spaces_ind:
            #        if i > ind_magmom_start:
            #            ind_magmom_end = i
            #            break
            #    magmom = float(line[ind_magmom_start:ind_magmom_end].split("=")[1])
            #    dict_energies["magmom_" + str(ind)] = magmom

        else:
            # print(line_split    )
            force_comps = np.array([float(i) for i in line_split[4:7]])
            # print(force_comps)
            forces.append(np.linalg.norm(force_comps))

    dict_energies["force"].append(np.mean(forces))  # gets last force instance
    return dict_energies


def main():
    directories = [
        # "./dft_full_slab_data/",
        # "./dft_full_slab_data/gas_phase/",
        # "./dft_full_slab_data/MD_NVE/",
        # "./dft_full_slab_data/fully_optimized_struc/"
        "/home/santiagovargas/dev/ml_PES_utils/will_data_061723/train/"
    ]

    # crawl all files in directory for files ending with .xyz
    for dir in directories:
        for root, dirs, files in os.walk(dir, topdown=True):
            for file in files:
                if file.endswith(".xyz"):
                    # check if csv exists
                    if not os.path.exists(root + "/" + file[:-4] + ".csv"):
                        print("Found xyz file: ", file)
                        lines = read_xyz(file=os.path.join(root, file))
                        dict_energies = parse_lines(lines)

                        df = pd.DataFrame(
                            {
                                "energy": dict_energies["energy"],
                                "force": dict_energies["force"],
                            }
                        )

                        energies = df["energy"].values
                        forces = df["force"].values
                        # get 75 percentile of energies and forces
                        energy_75 = np.percentile(energies, 75)
                        forces_75 = np.percentile(forces, 75)
                        energy_25 = np.percentile(energies, 25)
                        forces_25 = np.percentile(forces, 25)
                        iqr_energy = energy_75 - energy_25
                        iqr_forces = forces_75 - forces_25
                        iqr_range = 1.5
                        energy_outliers = list(
                            (energies < energy_25 - iqr_range * iqr_energy)
                            | (energies > energy_75 + iqr_range * iqr_energy)
                        )
                        if "bh" in file:
                            print("bh outlier finding")
                            force_outliers = (
                                list(
                                    (forces < forces_25 - iqr_range * iqr_forces)
                                    | (forces > forces_75 + iqr_range * iqr_forces)
                                )
                                | (forces > 3)
                                | (forces < -3)
                            )
                            # randomly select values not in the outliers to be outliers
                            for i in range(0, 12000):
                                rand_ind = np.random.randint(0, len(force_outliers))
                                force_outliers[rand_ind] = True
                        else:
                            force_outliers = list(
                                (forces < forces_25 - iqr_range * iqr_forces)
                                | (forces > forces_75 + iqr_range * iqr_forces)
                            )

                        df["force_outlier"] = force_outliers
                        df["energy_outlier"] = energy_outliers
                        print(file[:-4])
                        df.to_csv("./outlier_stats/" + file[:-4] + ".csv")


main()
