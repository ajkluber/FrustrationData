import os
import numpy as np
import matplotlib.pyplot as plt

import project_plotter
import project_util


class Dataset(project_util.Dataset):
    def __init__(self, *args, **kwargs):
        project_util.Dataset.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        data = np.load(self.datapath)
        with open("msm/Neff.dat", "w") as fout:
            fout.write("{}".format(data[0]))
        return data

if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]
    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]
    float_b = [ float(x) for x in b_values ]
    replicas = range(1,11)

    plotstyle = project_plotter.global_plot_style()

    datapath = "msm/Neff.npy"
    Neff_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$N_{eff}$", 
            "title":"$N_{eff}$ of MSM", "legend_loc":3,
            "saveas":"repavg_msm_Neff", 
            "newfig":True}
    Neff_plotspecs.update(plotstyle)

    Neff_data = Dataset(topologies, top_names, b_values, replicas, datapath)
    project_plotter.plot_data(Neff_data, Neff_plotspecs)

    plt.show()
