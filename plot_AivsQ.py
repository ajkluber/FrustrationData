import os
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter

class DatasetXvsY(project_util.DatasetXvsY):
    def __init__(self, *args, **kwargs):
        project_util.DatasetXvsY.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        data_x = np.loadtxt(self.datapath_x)
        data_y = np.load(self.datapath_y)
        #M = float(len(np.loadtxt("Ai_vs_Qtanh_0_05/nonnative_pairs.dat")))
        #data_y *= M
        return data_x, data_y
        

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
    datapath_x = "Ai_vs_" + coordname + "/Q_mid_bin.dat"
    datapath_y = "Ai_vs_" + coordname + "/A_vs_Q.npy"
    grid_dims = (4, 4)
    ylims = (0,0.07)

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_AvsQtanh_{}".format(top_names[i][j]))
        saveas.append(temp)

    A_plotspecs = {"xlabel":"Native Contacts", "ylabel":"A", 
            "title":"{} Non-native contacts", "saveas":saveas,
           
            "grid_dims":grid_dims, "coloridxs":coloridxs, "ylims":ylims, "xlims":(0,1),
            "xytext":(0.6,0.85)}

    A_grid_plotspecs = {"xlabel":"Native Contacts $Q$", "ylabel":"A", 
            "title":"Non-native contacts vs Q", "saveas":"repavg_grid_AvsQtanh",
           
            "grid_dims":(2,3), "coloridxs":coloridxs, "ylims":ylims, "xlims":(0,1),
            "xytext":(0.1,0.85)}

    A_raw_data = DatasetXvsY(topologies, top_names, b_values, replicas, datapath_x,
                datapath_y)

    A_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                for rep in range(len(A_raw_data.ydata[t][n][b])):
                    A_data.ydata[t][n][b].append(A_raw_data.ydata[t][n][b][rep])
                    A_data.xdata[t][n][b].append(A_raw_data.xdata[t][n][b][rep]/float(prot_n_native[t][n]))
    A_data._calc_repavg()

    project_plotter.plot_repavg(A_data, A_plotspecs)
    project_plotter.plot_repavg_protein_grid(A_data, A_grid_plotspecs)

    plt.show()
