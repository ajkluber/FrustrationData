import os
import glob
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter

class DatasetXvsY(project_util.DatasetXvsY):
    def __init__(self, *args, **kwargs):
        self.unfolded_Q = None
        self.n_res = None
        self.non_pairs = None
        project_util.DatasetXvsY.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        if self.unfolded_Q is None:
            # determine which Q bin to use
            U = np.loadtxt("Qtanh_0_05_profile/minima.dat")[0]
            self.unfolded_Q = U
        else:
            U = self.unfolded_Q

        if self.n_res is None:
            with open(glob.glob("T_*/ref.pdb")[0], "r") as fin:
                self.n_res = len(fin.readlines()) - 1
        if self.non_pairs is None:
            self.non_pairs = np.loadtxt("Ai_vs_Qtanh_0_05/nonnative_pairs.dat", dtype=int)

        Q_mid_bin = np.loadtxt("Ai_vs_Qtanh_0_05/Q_mid_bin.dat")
        Q_bin_idx = np.argmin((Q_mid_bin - U)**2)

        Ai_vs_Q = np.load("Ai_vs_Qtanh_0_05/Ai_vs_Q.npy")
        C = np.nan*np.zeros((self.n_res, self.n_res))
        for i in range(Ai_vs_Q.shape[1]):
            C[self.non_pairs[i,1], self.non_pairs[i,0]] = Ai_vs_Q[Q_bin_idx, i]
        # average over residue
        C_ma = np.ma.array(C, mask=np.isnan(C))

        data_x = np.arange(self.n_res)
        data_y = np.mean(C_ma, axis=1)
        return data_x, data_y

def get_Ai_contact_map():
    import numpy as np
    import matplotlib.pyplot as plt 

    #Q_mid_bin = np.loadtxt("Q_mid_bin.dat")

    # contact map
    C = np.nan*np.zeros((n_res, n_res))
    for i in range(n_non_dim):
        C[non_pairs[i,1], non_pairs[i,0]] = Ai_vs_Q[n, i]
    C_ma = np.ma.array(C, mask=np.isnan(C))

    plt.pcolormesh()

    return C_ma

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
    datapath_x = "Ai_vs_" + coordname + "/loop_length.dat"
    datapath_y = "Ai_vs_" + coordname + "/A_vs_loop_vs_Q.npy"
    grid_dims = (4, 4)
    #ylims = [[(0,0.1)],[(0,0.1)],[(0,0.1)],[(0,0.1)],[(0,0.1)]]
    ylims = (2e-4,1e-1)
    ylog = False
    xytext = (0.1, 0.9) 

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_AvsIndex_{}".format(top_names[i][j]))
        saveas.append(temp)

    plotspecs = {"xlabel":"Residue index", "ylabel":"A", "title":"{} Non-native contacts",
            "saveas":saveas,
            "grid_dims":grid_dims,
            "coloridxs":coloridxs, "xytext":xytext, "ylims":ylims,
            "ylog":ylog}

    dataset = DatasetXvsY(topologies, top_names, b_values, 
                replicas, datapath_x, datapath_y)

    project_plotter.plot_repavg(dataset, plotspecs)

    plt.show()
