import os
import numpy as np
import matplotlib.pyplot as plt

import project_util, project_plotter

if __name__ == "__main__":
    #topologies = ["alpha"]
    #top_names = [["1imq"]]

    #topologies = ["beta", "alpha"]
    #top_names = [["2akk"], ["1imq"]]

    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    replicas = range(1,11)

    plotstyle = project_plotter.global_plot_style()

    Tf_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"Temperature",
            "legend_loc":6, "saveas":None}
    Tf_plotspecs.update(plotstyle)

    Tg_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"Temperature",
            "saveas":"Tf_and_Tg"}
    Tg_plotspecs.update(plotstyle)
    Tg_plotspecs["newfig"] = False
    Tg_plotspecs["linestyle"] = "--"

    Tg_over_Tf_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":r"$\frac{T_g}{T_f}$", 
            "legend_loc":2, "saveas":"Tg_over_Tf"}
    Tg_over_Tf_plotspecs.update(plotstyle)

    Tf_data = project_util.Dataset(topologies, top_names, b_values, replicas, "Qtanh_0_05_profile/T_used.dat")
    Tg_data = project_util.Dataset(topologies, top_names, b_values, replicas, "Tg_calc/Tg_Enonnative.dat")

    Tg_over_Tf_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                Tf_b = Tf_data.data[t][n][b]
                Tg_b = Tg_data.data[t][n][b]
                for rep in range(len(Tf_b)):
                    if not np.isnan(Tf_b[rep]) and not np.isnan(Tg_b[rep]):
                        Tg_over_Tf_data.data[t][n][b][rep] = Tg_b[rep]/Tf_b[rep]
    Tg_over_Tf_data._calc_repavg()

    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            temp_std = []
            for b in range(len(b_values)):
                Tf_b = Tf_data.data[t][n][b]
                Tg_b = Tg_data.data[t][n][b]
                good1 = ~np.isnan(Tg_b)
                good2 = ~np.isnan(Tf_b)
                if (np.sum(good1) > 1) and (np.sum(good2) > 1):
                    da_over_a = (np.std(Tg_b[good1])/np.sqrt(len(Tg_b[good1])))/(np.mean(Tg_b[good1]))
                    dx_over_x = (np.std(Tf_b[good2])/np.sqrt(len(Tf_b[good2])))/(np.mean(Tf_b[good2]))
                    Q = np.mean(Tg_b[good1])/np.mean(Tf_b[good2])
                    dQ = Q*np.sqrt(da_over_a**2 + dx_over_x**2)
                    temp_std.append(dQ)
                else:
                    temp_std.append(np.nan)
            Tg_over_Tf_data.stddata[t][n] = np.array(temp_std)

    project_plotter.plot_data(Tf_data, Tf_plotspecs, err_mult=1.96)
    project_plotter.plot_data(Tg_data, Tg_plotspecs, err_mult=1.96)
    project_plotter.plot_data(Tg_over_Tf_data, Tg_over_Tf_plotspecs, err_mult=1.96)
    plt.show()
