# use csv to merge trajectories
import os
import pandas as pd


def read_xyz(file="./gas_phase/C3H6/c3h6.xyz"):
    """
    Read the xyz file and return the coordinates of the atoms
    """
    with open(file, "r") as f:
        lines = f.readlines()
    return lines


def parse_lines(lines, write_to, outlier_force, outlier_energy):
    ind = -1
    outlier_count = 0
    # start_control = False
    for i, line in enumerate(lines):
        line_split = line.split()
        if line_split[0].isdigit():
            if i < len(lines) - 1:
                # if (lines[i+1].split("=")[0] == "Properties"):
                #    properties_tag = True
                #
                # else:
                #    properties_tag = False
                ind += 1
                if outlier_force[ind] == 1 or outlier_energy[ind] == 1:
                    outlier_count += 1
        if outlier_force[ind] == False and outlier_energy[ind] == False:
            write_to.write(line)

    print("Outliers: ", outlier_count)
    print("Total: ", ind)


def main():
    # get all files in directory
    csv = os.listdir()

    directories = [
        # "./dft_full_slab_data/",
        # "./dft_full_slab_data/gas_phase/",
        # "./dft_full_slab_data/MD_NVE/",
        # "./dft_full_slab_data/fully_optimized_struc/"
        # "./dft_full_slab_data/transf_set_al2o3",
        "/home/santiagovargas/dev/ml_PES_utils/will_data_061723/train/"
        # "./dft_full_slab_data/wrong_data/"
    ]

    # open master xyz file
    master_xyz = open("target_md_bh_aug_full.xyz", "w")

    for dir in directories:
        for root, dirs, files in os.walk(dir, topdown=True):
            for file in files:
                if file.endswith(".xyz"):
                    # check if csv exists
                    #print("./outlier_stats/" + file[:-4] + ".csv")
                    #if os.path.exists("./outlier_stats/" + file[:-4] + ".csv"):
                        print("using csv to merge trajectories")
                        lines = read_xyz(file=os.path.join(root, file))
                        df = pd.read_csv("./outlier_stats/" + file[:-4] + ".csv")
                        force_outliers = df["force_outlier"].values
                        energy_outliers = df["energy_outlier"].values

                        # if not an outlier, append line to list
                        parsed_lines = parse_lines(
                            lines, master_xyz, force_outliers, energy_outliers
                        )

                    else:
                        print("No csv file found, skipping")


main()
