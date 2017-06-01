import os
import glob
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter

class DatasetXvsY(project_util.DatasetXvsY):
    def __init__(self, *args, **kwargs):
        project_util.DatasetXvsY.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        data_y = np.abs(np.load(self.datapath_y))
        data_x = np.arange(1, len(data_y) + 1)
        return data_x, data_y

if __name__ == "__main__":
    topologies = ["alpha","beta", "mixed"]
    top_names = [["1r69", "1imq"],["1fmk", "2akk"], ["1e0g"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]
    coloridxs = [ b/1.3 for b in float_b ]

    replicas = range(1,11)

    coordname = "Qtanh_0_05"
    datapath_x = "persist_length/corr_ij.npy"
    datapath_y = "persist_length/corr_ij.npy"
    grid_dims = (4, 4)
    #ylims = [[(0,0.1)],[(0,0.1)],[(0,0.1)],[(0,0.1)],[(0,0.1)]]
    ylims = (1e-4, 1)
    ylog = True
    xytext = (0.6, 0.9) 

    saveas = []
    for i in range(len(topologies)):
        temp = []
        for j in range(len(top_names[i])):
            temp.append("repavg_corrij_{}".format(top_names[i][j]))
        saveas.append(temp)

    plotspecs = {"xlabel":"Separation $l$", 
            "ylabel":"$\\langle r_i \cdot r_{i+ l}\\rangle$", "title":"{} correlation",
            "saveas":saveas,
            "grid_dims":grid_dims,
            "coloridxs":coloridxs, "xytext":xytext, "ylims":ylims,
            "ylog":ylog}

    dataset = DatasetXvsY(topologies, top_names, b_values, 
                replicas, datapath_x, datapath_y)

    project_plotter.plot_repavg(dataset, plotspecs)

    plt.show()
