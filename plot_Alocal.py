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
        Q_mid_bin = np.loadtxt("Alocal_vs_Qtanh_0_05/Q_mid_bin.dat")
        Q_bin_idx = np.argmin((Q_mid_bin - U)**2)


        data_y = np.load(self.datapath_y)[Q_bin_idx,:]
        data_x = np.loadtxt(self.datapath_x)
        #data_y = np.load(self.datapath_y)
        #M = float(len(np.loadtxt("Alocal_vs_Qtanh_0_05/nonnative_pairs.dat")))
        #data_y *= M
        return data_x, data_y

if __name__ == "__main__":
    topologies = ["alpha","beta", "mixed"]
    top_names = [["1r69", "1imq"],["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    prot_n_native = [[155, 219], [164, 213], [117]]
    #N = np.array([63, 86, 58, 74, 48])
    #M = np.array([155, 219, 164, 213, 117])

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]
    coloridxs = [ b/1.3 for b in float_b ]

    replicas = range(1,11)

    coordname = "Qtanh_0_05"
    datapath_x = "Alocal_vs_" + coordname + "/Q_mid_bin.dat"
    datapath_y = "Alocal_vs_" + coordname + "/Alocal_vs_Q.npy"
    grid_dims = (2, 3)
    xytext = (0.1, 0.8) 

    Alocal_plotspecs = {"xlabel":"Q", "ylabel":"$A_{local}$ / N", 
            "title":"Non-native contacts",
            "saveas":"repavg_grid_Alocal_vs_Qtanh",
            "grid_dims":grid_dims, "ylims":(0,60./48),
            "coloridxs":coloridxs, "xytext":xytext}

    Alocal_data = project_util.DatasetXvsY(topologies, top_names, b_values, 
                replicas, datapath_x, datapath_y)

    # plot maximum A scaled by M, the number of native contacts.
    Alocal2_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                for rep in range(len(Alocal_data.ydata[t][n][b])):
                    Alocal2_data.ydata[t][n][b].append(Alocal_data.ydata[t][n][b][rep]/float(prot_sizes[t][n]))
                    Alocal2_data.xdata[t][n][b].append(Alocal_data.xdata[t][n][b][rep]/float(prot_n_native[t][n]))
    Alocal2_data._calc_repavg()

    Anonlocal_plotspecs = {"xlabel":"Q", "ylabel":"$A_{non-local} / N$", 
            "title":"Non-native contacts",
            "saveas":"repavg_grid_Anonlocal_vs_Qtanh",
            "grid_dims":grid_dims, "ylims":(0,60./48),
            "coloridxs":coloridxs, "xytext":xytext}

    Anonlocal_data = project_util.DatasetXvsY(topologies, top_names, b_values, 
                replicas, datapath_x, "Alocal_vs_Qtanh_0_05/Anonlocal_vs_Q.npy")

    # plot maximum A scaled by M, the number of native contacts.
    Anonlocal2_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                for rep in range(len(Anonlocal_data.ydata[t][n][b])):
                    Anonlocal2_data.ydata[t][n][b].append(Anonlocal_data.ydata[t][n][b][rep]/float(prot_sizes[t][n]))
                    Anonlocal2_data.xdata[t][n][b].append(Anonlocal_data.xdata[t][n][b][rep]/float(prot_n_native[t][n]))
    Anonlocal2_data._calc_repavg()

    project_plotter.plot_repavg_protein_grid(Alocal2_data, Alocal_plotspecs)
    project_plotter.plot_repavg_protein_grid(Anonlocal2_data, Anonlocal_plotspecs)

    plt.show()
