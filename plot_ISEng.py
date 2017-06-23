import os
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter

import inherent_structure.util as util

class DatasetXvsY(project_util.DatasetXvsY):
    def __init__(self, *args, **kwargs):
        self._get_data = kwargs.pop("data_getter")
        project_util.DatasetXvsY.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        data_x, data_y = self._get_data()
        return data_x, data_y

def get_Enon_hist():
    nbins_Enat = 70
    nbins_Enon = 20
    Enat, Enon = util.get_native_nonnative_energies()
    P_Enat, Enat_mid_bin, U = util.determine_U_frames(Enat, nbins_Enat)

    P_Enon_U, Enon_mid_bin = np.histogram(Enon[U], bins=nbins_Enon)
    return Enon_mid_bin, P_Enon_U


if __name__ == "__main__":
    topologies = ["alpha","beta", "mixed"]
    top_names = [["1r69", "1imq"],["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    prot_n_native = [[155, 219], [164, 213], [117]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25",
                "1.35", "1.45"]

    float_b = [ float(x) for x in b_values ]
    coloridxs = [ b/1.3 for b in float_b ]

    replicas = range(1,11)

    coordname = "Qtanh_0_05"
    datapath_x = "Qtanh_0_05_profile/T_used.dat"
    datapath_y = "T_{}_1/inherent_structures/Enat.npy"
    grid_dims = (4, 4)
    ylims = (0,0.07)


    Enon_plotspecs = {"xlabel":"Native Contacts", "ylabel":"A", 
            "title":"{} Non-native contacts", "saveas":None,
            "grid_dims":grid_dims, "ylims":ylims, "xlims":(0,1),
            "xytext":(0.6,0.85)}

    Enon_grid_plotspecs = {"xlabel":"Native Contacts $Q$", "ylabel":r"$E_{non}$",
            "title":"Non-native contacts vs Q", "saveas":None,
            "grid_dims":(2,3), "coloridxs":coloridxs, "ylims":ylims, "xlims":(0,1),
            "xytext":(0.1,0.85)}

    P_Enon_data = DatasetXvsY(topologies, top_names, b_values, replicas, datapath_x,
                datapath_y, T_used=True, data_getter=get_Enon_hist)

