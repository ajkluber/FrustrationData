import os
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter


class Dataset(project_util.Dataset):

    def __init__(self, *args, **kwargs):
        self.REM_idx = kwargs.pop("REM_idx")
        project_util.Dataset.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        with open(self.datapath, "r") as fin:
            data = float(fin.read().split()[self.REM_idx])
        return data

if __name__ == "__main__":
    #topologies = ["alpha"]
    #top_names = [["1imq"]]

    #topologies = ["beta", "alpha"]
    #top_names = [["2akk"], ["1imq"]]

    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    markers = ["^", "o", "s"]

    legend_key = []
    for t in range(len(topologies)):
        top_keys = []
        for n in range(len(top_names[t])):
           top_keys.append("{} ({})".format(top_names[t][n], prot_sizes[t][n]))
        legend_key.append(top_keys)
    #ylims = (0,1)

    color_codes = [["#5DA5DA", "#F17CB0"], ["#60BD68", "#B2912F"], ["#FAA43A"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    replicas = range(1,11)

    dEstab_data = Dataset(topologies, top_names, b_values, replicas, "Tg_calc/REM_parameters.dat", REM_idx=0)
    dEnon_data = Dataset(topologies, top_names, b_values, replicas, "Tg_calc/REM_parameters.dat", REM_idx=1)
#    Ebar_data = Dataset(topologies, top_names, b_values, replicas, "Tg_calc/REM_parameters.dat", REM_idx=2)
#    S0_data = Dataset(topologies, top_names, b_values, replicas, "Tg_calc/REM_parameters.dat", REM_idx=3)

    dEstab_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$\\delta E_{stab}$ ($k_BT$)",
            "title":"Stability gap", "legend_loc":3, "ylims":(0,160),
            "saveas":"dEstab_vs_b", 
           }

    dEnon_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$\\Delta E_{non}$ ($k_BT$)",
            "title":"Roughness", "legend_loc":2, "ylims":(0,18),
            "saveas":"dEnon_vs_b", 
           }

    Ebar_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$\\overline{E_{non}}$ ($k_BT$)",
            "title":"Mean energy", "legend_loc":2, "ylims":(0,230),
            "saveas":"Ebar_vs_b", 
           }

    S0_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$S_0$ ($k_B$)",
            "title":"Entropy", "legend_loc":3, "ylims":(0,150),
            "saveas":"S0_vs_b", 
           }

#    # Plot Random Energy Model parameters 
#    project_plotter.plot_data(dEstab_data, dEstab_plotspecs)
#    project_plotter.plot_data(Ebar_data, Ebar_plotspecs)
#    project_plotter.plot_data(S0_data, S0_plotspecs)
#    project_plotter.plot_data(dEnon_data, dEnon_plotspecs)
#    #plt.ylim(0,18**2)
#    #b = np.linspace(0, 1.4, 1000)
#    #plt.plot(b, 11*b*b, 'k', lw=2)
#    #plt.plot(b, 145*(b**4), 'k', lw=2)
#    #plt.plot(b, 50*(b**4), 'k', lw=2)
#    #plt.plot(b, 150*(b**3), 'k', lw=2)
#
#    # plot roughness scaled by size
#    dEnon_scaled_data = project_util.Dataset(topologies, top_names, b_values, replicas,"", empty=True)
#    for t in range(len(topologies)):
#        for n in range(len(top_names[t])):
#            for b in range(len(b_values)):
#                dEnon_data_for_b = dEnon_data.data[t][n][b]
#                for rep in range(len(dEnon_data_for_b)):
#                    dEnon_scaled_data.data[t][n][b].append(dEnon_data_for_b[rep]/float(prot_sizes[t][n]))
#
#    dEnon_scaled_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$\\Delta E_{non}$ / N (kJ/mol)",
#            "title":"Roughness", "legend_loc":2,
#            "saveas":"dEnon_scaled_vs_b", 
#           }
#
#    project_plotter.plot_data(dEnon_scaled_data, dEnon_scaled_plotspecs)

    # plot roughness scaled by size
    dEnon_over_dEstab_data = project_util.Dataset(topologies, top_names, b_values, replicas,"", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                dEnon = dEnon_data.data[t][n][b]
                dEstab = dEstab_data.data[t][n][b]
                for rep in range(len(dEnon)):
                    dEnon_over_dEstab_data.data[t][n][b].append(dEnon[rep]/dEstab[rep])

    dEnon_over_dEstab_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$\\Delta E_{non}$ / $\\delta E_{stab}$",
            "title":" ", "legend_loc":2, "ylims":(0,0.5),
            "saveas":"dEnon_over_dEstab_vs_b", 
           }

    project_plotter.plot_data(dEnon_over_dEstab_data, dEnon_over_dEstab_plotspecs)


    plt.show()
