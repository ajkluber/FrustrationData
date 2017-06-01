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
    datapath_x = coordname + "_profile/T_{:.2f}_mid_bin.dat"
    datapath_y = coordname + "_profile/T_{:.2f}_F.dat"
    grid_dims = (4, 4)
    ylims = (0., 5)
    xytext = (0.1, 0.9) 

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_FvsQtanh_{}".format(top_names[i][j]))
        saveas.append(temp)

    F_plotspecs = {"xlabel":"Native Contacts", "ylabel":"Free Energy (k$_B$T)",
            "title":"{} Free energy profiles", "saveas":saveas,
            "grid_dims":grid_dims, "ylims":ylims, "coloridxs":coloridxs}

    F_grid_plotspecs = {"xlabel":"Native contacts Q", "ylabel":"Free Energy (k$_B$T)", 
            "title":"Free energy profiles",
            "saveas":"repavg_grid_FvsQtanh",
            "grid_dims":(2,3),
            "coloridxs":coloridxs, "xytext":xytext, "ylims":(0, 6), "xlims":(0,1)}

    F_raw_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas, datapath_x,
            datapath_y, T_used=True)

    F_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                for rep in range(len(F_raw_data.ydata[t][n][b])):
                    F_data.xdata[t][n][b].append(F_raw_data.xdata[t][n][b][rep]/float(prot_n_native[t][n]))
                    F_data.ydata[t][n][b].append(F_raw_data.ydata[t][n][b][rep])
    F_data._calc_repavg()

    project_plotter.plot_repavg(F_data, F_plotspecs)

    project_plotter.plot_repavg_protein_grid(F_data, F_grid_plotspecs)

    plt.show()
