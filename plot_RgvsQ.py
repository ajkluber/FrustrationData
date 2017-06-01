import os
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter

if __name__ == "__main__":
    topologies = ["alpha","beta", "mixed"]
    top_names = [["1r69", "1imq"],["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    prot_n_native = [[155, 219], [164, 213], [117]]

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
    xytext = (0.1, 0.9) 

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_RgvsQtanh_{}".format(top_names[i][j]))
        saveas.append(temp)

    Rg_plotspecs = {"xlabel":"Native contacts Q", "ylabel":"$R_g$ (nm)", "title":"{} Radius of gyration",
            "saveas":saveas,
            "grid_dims":grid_dims,
            "coloridxs":coloridxs, "xytext":xytext, "ylims":ylims}

    Rg_grid_plotspecs = {"xlabel":"Native contacts Q", "ylabel":"$R_g$ (nm)", 
            "title":"Radius of gyration",
            "saveas":"repavg_grid_RgvsQtanh",
            "grid_dims":(2,3),
            "coloridxs":coloridxs, "xytext":xytext, "ylims":ylims, "xlims":(0,1)}

    Rg_raw_data = project_util.DatasetXvsY(topologies, top_names, b_values, 
                replicas, datapath_x, datapath_y)

    Rg_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                for rep in range(len(Rg_raw_data.ydata[t][n][b])):
                    Rg_data.xdata[t][n][b].append(Rg_raw_data.xdata[t][n][b][rep]/float(prot_n_native[t][n]))
                    Rg_data.ydata[t][n][b].append(Rg_raw_data.ydata[t][n][b][rep])
    Rg_data._calc_repavg()

    #project_plotter.plot_repavg(Rg_data, Rg_plotspecs)

    project_plotter.plot_repavg_protein_grid(Rg_data, Rg_grid_plotspecs)



    plt.show()
