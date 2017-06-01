import os
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter

class DatasetXvsY(project_util.DatasetXvsY):
    def __init__(self, *args, **kwargs):
        self.unfolded_Q = None
        project_util.DatasetXvsY.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        if self.unfolded_Q is None:
            # determine which Q bin to use
            U = np.loadtxt("Qtanh_0_05_profile/minima.dat")[0]
            self.unfolded_Q = U
        else:
            U = self.unfolded_Q
        U = np.loadtxt("Qtanh_0_05_profile/minima.dat")[0]
        Q_mid_bin = np.loadtxt("Ai_vs_Qtanh_0_05/Q_mid_bin.dat")
        Q_bin_idx = np.argmin((Q_mid_bin - U)**2)
        data_y = np.load(self.datapath_y)[Q_bin_idx,:]
        data_x = np.loadtxt(self.datapath_x)
        #data_y = np.load(self.datapath_y)
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
    datapath_x = "Ai_vs_" + coordname + "/loop_length.dat"
    datapath_y = "Ai_vs_" + coordname + "/A_vs_loop_vs_Q.npy"
    grid_dims = (4, 4)
    #ylims = (2e-4,1e-1)
    #ylog = False
    #xytext = (0.1, 0.9) 
    ylims = (2e-4,0.07)
    avg_ylims = (0, 0.063)
    ylog = False
    xytext = (0.2, 0.85) 

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_AvsLoop_{}".format(top_names[i][j]))
        saveas.append(temp)

    AvsLoop_plotspecs = {"xlabel":"Sequence separation", "ylabel":"$A$", 
            "title":"{} Non-native contacts",
            "saveas":saveas,
            "grid_dims":grid_dims,
            "coloridxs":coloridxs, "xytext":xytext, "ylims":ylims,
            "ylog":ylog, "avg_ylims":avg_ylims}

    AvsLoop_grid_plotspecs = {"xlabel":"Sequence separation", "ylabel":"$A$", 
            "title":" ",
            "saveas":"repavg_grid_AvsLoop",
            "grid_dims":(2,3),
            "coloridxs":coloridxs, "xytext":xytext, "ylims":(0, 0.1)}

    AvsLoop_data = DatasetXvsY(topologies, top_names, b_values, 
                replicas, datapath_x, datapath_y)

    project_plotter.plot_repavg(AvsLoop_data, AvsLoop_plotspecs)
    project_plotter.plot_repavg_protein_grid(AvsLoop_data, AvsLoop_grid_plotspecs)

    plt.show()
