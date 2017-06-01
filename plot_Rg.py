import os
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter

if __name__ == "__main__":
    topologies = ["alpha","beta", "mixed"]
    top_names = [["1r69", "1imq"],["1fmk", "2akk"], ["1e0g"]]
    #topologies = ["alpha", "beta"]
    #top_names = [["1imq"], ["2akk"]]
    #topologies = ["alpha"]
    #top_names = [["1r69"]]
    #topologies = ["beta"]
    #top_names = [["1fmk"]]
    #topologies = ["mixed"]
    #top_names = [["1e0g"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]
    coloridxs = [ b/1.3 for b in float_b ]

    replicas = range(1,11)

    coordname = "Qtanh_0_05"
    datapath_x = "Rg_vs_" + coordname + "/mid_bin.npy"
    datapath_y = "Rg_vs_" + coordname + "/Rg_vs_bin.npy"
    grid_dims = (4, 4)
    ylims = (0.5, 2.5)
    ylog = False
    xytext = (0.1, 0.9) 

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_RgvsQtanh_{}".format(top_names[i][j]))
        saveas.append(temp)

    plotspecs = {"xlabel":"Native contacts Q", "ylabel":"$R_g$ (nm)", "title":"{} Radius of gyration",
            "saveas":saveas,
            "grid_dims":grid_dims,
            "coloridxs":coloridxs, "xytext":xytext, "ylims":ylims,
            "ylog":ylog}

    dataset = project_util.DatasetXvsY(topologies, top_names, b_values, 
                replicas, datapath_x, datapath_y)

    project_plotter.plot_repavg(dataset, plotspecs)

    plt.show()
