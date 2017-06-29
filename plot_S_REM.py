import os
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter


class Dataset(project_util.Dataset):

    def __init__(self, *args, **kwargs):
        self.REM_idx = kwargs.pop("REM_idx")
        project_util.Dataset.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        with open(self.datapath, "r") as fin:
            data = float(fin.read().split()[self.REM_idx])
        return data

if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25",
                "1.35", "1.45"]

    float_b = [ float(x) for x in b_values ]
    coloridxs = [ b/1.3 for b in float_b ]

    prot_sizes = [[63, 86], [58, 74], [48]]

    replicas = range(1,11)

    # get data
    S_REM_data = project_util.DatasetXvsY(topologies, top_names, 
            b_values, replicas, "Tg_calc/Enon_mid_bin.npy", "Tg_calc/S_REM_Enon.npy")

    plotstyle = project_plotter.global_plot_style()

    grid_dims = (4, 4)
    ylims = (-20, 200)
    xytext = (0.1, 0.9) 

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_S_REM_vs_Enon_{}".format(top_names[i][j]))
        saveas.append(temp)

    S_REM_plotspecs = {"xlabel":"Native Contacts", 
            "ylabel":r"$S_{REM}(E_{non})$ ($k_B$)",
            "title":"{} Entropy", "saveas":None, "ylims":ylims,
            "grid_dims":grid_dims, "coloridxs":coloridxs}
    S_REM_plotspecs.update(plotstyle)

    #S_REM_plotspecs = {"xlabel":"Frustration ($b$)", 
    #        "ylabel":r"$S_0$ ($k_B$)",
    #        "title":"Entropy", "legend_loc":3, 
    #        "saveas":"S0_vs_b"}
    #S_REM_plotspecs.update(plotstyle)

    project_plotter.plot_repavg(S_REM_data, S_REM_plotspecs)
    plt.show()
