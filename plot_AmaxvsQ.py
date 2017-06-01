import os
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter

if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"],["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    prot_n_native = [[155, 219], [164, 213], [117]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]
    coloridxs = [ b/1.3 for b in float_b ]

    replicas = range(1, 11)

    coordname = "Qtanh_0_05"
    datapath_x = "Amax_vs_" + coordname + "/Q_mid_bin.dat"
    datapath_y = "Amax_vs_" + coordname + "/Amax_vs_Q.npy"
    grid_dims = (4, 4)
    #ylims = [[(0, 0.07)], [(0., 0.07)], [(0., 0.07)]]
    #ylims = (0, 0.07)

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_AmaxvsQtanh_{}".format(top_names[i][j]))
        saveas.append(temp)

    Amax_plotspecs = {"xlabel":"Native Contacts", "ylabel":"$A_{max}$", 
            "title":"{} Non-native contacts", "saveas":saveas,
           
            "grid_dims":grid_dims, "coloridxs":coloridxs, 
            "xytext":(0.6,0.85)}

    Amax_grid_plotspecs = {"xlabel":"Native Contacts", "ylabel":"$A_{max}$", 
            "title":" ", "saveas":"repavg_grid_Amax_vs_Qtanh",
           
            "grid_dims":(2,3), "coloridxs":coloridxs, 
            "xytext":(0.6,0.85)}

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_Amax_scaled_N_vsQtanh_{}".format(top_names[i][j]))
        saveas.append(temp)

    Amax_scaled_N_plotspecs = {"xlabel":"Native Contacts", "ylabel":"$A_{max} / N$", 
            "title":"{} Non-native contacts", "saveas":saveas, "ylims":(0, 2.5), "avg_ylims":(0,2.5),
           
            "grid_dims":grid_dims, "coloridxs":coloridxs, 
            "xytext":(0.6,0.85)}

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_Amax_scaled_M_vsQtanh_{}".format(top_names[i][j]))
        saveas.append(temp)

    Amax_scaled_M_plotspecs = {"xlabel":"Native Contacts", "ylabel":"$A_{max} / M$", 
            "title":"{} Non-native contacts", "saveas":saveas, "ylims":(0, 1.1), "avg_ylims":(0, 1.1),
            "xlims":(0,1), "avg_xlims":(0,1), 
            "grid_dims":grid_dims, "coloridxs":coloridxs, 
            "xytext":(0.6,0.85)}

    # Get raw Amax data
    Amax_scaled_M_grid_plotspecs = {"xlabel":"Native Contacts", "ylabel":"$A_{max} / M$", 
            "title":" ", "saveas":"repavg_grid_Amax_scaled_M_vs_Qtanh", "ylims":(0, 1.1),
            "xlims":(0,1),
            "grid_dims":(2,3), "coloridxs":coloridxs, 
            "xytext":(0.6,0.85)}

    Amax_raw_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas, datapath_x, datapath_y)

    # plot maximum A scaled by N, the number of residues.
    Amax_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                for rep in range(len(Amax_raw_data.ydata[t][n][b])):
                    Amax_data.ydata[t][n][b].append(Amax_raw_data.ydata[t][n][b][rep])
                    Amax_data.xdata[t][n][b].append(Amax_raw_data.xdata[t][n][b][rep]/float(prot_n_native[t][n]))
    Amax_data._calc_repavg()

    # plot maximum A scaled by N, the number of residues.
    Amax_scaled_N_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                for rep in range(len(Amax_data.ydata[t][n][b])):
                    Amax_scaled_N_data.ydata[t][n][b].append(Amax_data.ydata[t][n][b][rep]/float(prot_sizes[t][n]))
                    Amax_scaled_N_data.xdata[t][n][b].append(Amax_data.xdata[t][n][b][rep])
    Amax_scaled_N_data._calc_repavg()



    # plot maximum A scaled by M, the number of native contacts.
    Amax_scaled_M_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                for rep in range(len(Amax_data.ydata[t][n][b])):
                    Amax_scaled_M_data.ydata[t][n][b].append(Amax_data.ydata[t][n][b][rep]/float(prot_n_native[t][n]))
                    Amax_scaled_M_data.xdata[t][n][b].append(Amax_data.xdata[t][n][b][rep])
    Amax_scaled_M_data._calc_repavg()


    #project_plotter.plot_repavg(Amax_data, Amax_plotspecs)
    #project_plotter.plot_repavg(Amax_scaled_N_data, Amax_scaled_N_plotspecs)
    #project_plotter.plot_repavg(Amax_scaled_M_data, Amax_scaled_M_plotspecs)

    project_plotter.plot_repavg_protein_grid(Amax_scaled_M_data, Amax_scaled_M_grid_plotspecs)
    project_plotter.plot_repavg_protein_grid(Amax_data, Amax_grid_plotspecs)

    plt.show()
