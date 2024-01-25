import os
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def read_xyz(file="./gas_phase/C3H6/c3h6.xyz"):
    """Read the xyz file and return the coordinates of the atoms"""
    with open(file, "r") as f:
        lines = f.readlines()
    return lines


def parse_lines(lines):
    # find lines that start with a number
    # and split them into a list of strings
    dict_energies = {"energy": [], "force": [], "stress": [], "lines": []}

    ind = -1
    start_control = False

    for ind, line in enumerate(lines):
        line_split = line.split()
        # print(line_split)

        if line_split[0].isdigit():
            lines_temp = [line]

        elif line[0:7] == "Lattice":
            lines_temp.append(line)
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

        else:
            lines_temp.append(line)
            force_comps = np.array([float(i) for i in line_split[4:7]])
            forces.append(np.linalg.norm(force_comps))

        if ind < len(lines) - 1:
            if lines[ind + 1].split()[0].isdigit():
                # if lines[ind+1][0:7] == "Lattice":
                dict_energies["lines"].append(lines_temp)
    dict_energies["lines"].append(lines_temp)
    dict_energies["force"].append(np.mean(forces))  # gets last force instance
    return dict_energies


def main():
    # file = './dft_full_slab_data/fully_optimized_struc/Al2O3+Pt8+Hx/al2o3_pt8_hx.xyz'
    file = "./target_md_bh_aug_full.xyz"
    # file = './merged_al2o3_10257.xyz'
    # file = './merged_tot_23693.xyz'
    # file = './dft_full_slab_data/perturbed/perturb.xyz'
    # file = './dft_full_slab_data/fully_optimized_struc/Al2O3_Pt8_Hx_C3H6/new_al2o3_pt8_hx_c3h6.xyz'
    # root = "./"

    lines = read_xyz(file=file)

    dict_energies = parse_lines(lines)

    train_test_split = 0.85
    # randomly draw train_test_split*100% of the data
    # and use the rest for testingcd
    print(len(dict_energies["energy"]))
    print(len(dict_energies["lines"]))

    train_ind = np.random.choice(
        len(dict_energies["energy"]),
        int(train_test_split * len(dict_energies["energy"])),
        replace=False,
    )
    test_ind = np.setdiff1d(np.arange(len(dict_energies["energy"])), train_ind)

    # test file as name of the file + _test.xyz
    # train file as name of the file + _train.xyz
    # print(file.split("."))
    test_file = "." + file.split(".")[1] + "_test.xyz"
    train_file = "." + file.split(".")[1] + "_train.xyz"
    # open file to write to - training
    with open(train_file, "w") as f:
        # write the lines to the file
        for i in train_ind:
            f.writelines(dict_energies["lines"][i])
    # open file to write to - testing
    with open(test_file, "w") as f:
        for i in test_ind:
            f.writelines(dict_energies["lines"][i])


main()
