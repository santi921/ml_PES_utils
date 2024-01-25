import os

# import ase
from ase import db
from glob import glob

# import matplotlib.pyplot as plt

import numpy as np
from ase_notebook import AseView, ViewConfig, concatenate_svgs

# from ase_notebook.backend.svg import svg_to_pdf
# import svgwrite


def plot_db_minima(
    database_file, max_plot=20, orientation_str="90x,90y,0z", figsize=(20, 100)
):
    print("plotting from {}".format(database_file))
    ase_view = AseView(
        rotations=orientation_str,
        atom_font_size=16,
        axes_length=30,
        canvas_size=(1000, 1000),
        zoom=1.2,
        show_bonds=True,
    )
    idx = -1
    count_plot = 0
    # create sub_plot
    ase_view.config.atom_show_label = False
    ase_view.config.canvas_background_opacity = 0.0
    svg_list = []
    energy_list = []
    # if max_plot > 0:
    #    rows = int(np.ceil(max_plot / cols))
    # else:
    #    rows = int(np.ceil(len(minima_dbs) / cols))

    # fig, axs = plt.subplots(rows, cols, figsize=figsize)
    # axs = axs.flatten()
    db_minima = db.connect(database_file)
    print("total structures: {}".format(db_minima.count()))
    # print columns of sqlite3 database
    if max_plot > 0 and max_plot < db_minima.count():
        # randomly select structure sto plot
        idx = np.random.choice(db_minima.count(), max_plot, replace=False)
        print("selected {} structures".format(len(idx)))
        print("type index {}".format(type(idx)))
    else:
        idx = np.arange(db_minima.count())

    db_sub = db_minima.select("age<2y")

    for i, row in enumerate(db_sub):
        atoms = row.toatoms()
        if type(idx) == np.ndarray:
            if i in idx:
                svg = ase_view.make_svg(atoms, center_in_uc=True)
                energy = row.energy
                # write energy on svg
                energy_list.append(energy)
                svg.add(
                    svg.text(
                        "{:.2f}".format(energy),
                        insert=(250, 150),
                        font_size="30px",
                        font_family="Arial",
                    )
                )
                # print(type(svg))
                svg_list.append(svg)
                count_plot += 1
        else:
            svg = ase_view.make_svg(atoms, center_in_uc=True)

            # write energy on svg
            energy = row.energy
            energy_list.append(energy)
            svg.add(
                svg.text(
                    "{:.2f}".format(energy),
                    insert=(250, 150),
                    font_size="30px",
                    font_family="Arial",
                )
            )
            svg_list.append(svg)

    # tight layout
    # sort by energy
    sort_idx = np.argsort(energy_list)
    svg_list = [svg_list[i] for i in sort_idx]
    print("SVGS: {}".format(len(svg_list)))
    return svg_list
    # return svg_list, concatenate_svgs(svg_list, max_columns=cols, scale=0.5, label=False)


def main():
    run = "2/low/"

    data_root = (
        "/home/santiagovargas/dev/ml_PES_utils/data/local_min_h_filter_1012/" + run
    )
    minima_dbs = glob(data_root + "*.db")

    for db in minima_dbs:
        raw_list = plot_db_minima(
            db,
            max_plot=250,
            orientation_str="270x,270y,0z",
            figsize=(20, 100),
        )

        save_folder = "../../data/figures_bh_filter_1012/" + run
        run_name = db.split("/")[-1].split(".")[0]

        print("model name: {}".format(run_name))
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)

        for i, svg in enumerate(raw_list):
            svg.saveas(save_folder + "model_{}_{}.svg".format(run_name, i))
            # svg_to_pdf(svg, save_folder + "model_{}_{}.pdf".format(run_name, i))
            print("saved model_{}_{}.svg".format(run_name, i))


main()
