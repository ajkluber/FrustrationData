import os
import numpy as np
import matplotlib.pyplot as plt

import project_plotter
import project_util

import plotter.pmfutil


if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]
    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25", "1.35", "1.45"]
    float_b = [ float(x) for x in b_values ]

    replicas = range(1, 11)

    plotstyle = project_plotter.global_plot_style()

    coordname = "Qtanh_0_05"
    datapath = coordname + "_transit_time/forward_mean"

    tp_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$t_{tp}$ (frames)", 
            "legend_loc":6,
            "saveas":"repavg_tp_1D", "saveas_formats":["png","pdf"]}
    tp_plotspecs.update(plotstyle)

    tp_norm_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$t_{tp}$ (frames)", 
            "legend_loc":4,
            "saveas":"repavg_tp_norm_1D", "saveas_formats":["png","pdf"]}
    tp_norm_plotspecs.update(plotstyle)

    tp_ylog_norm_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$t_{tp}$ (frames)", 
            "legend_loc":4,
            "saveas":"repavg_tp_ylog_norm_1D", "saveas_formats":["png","pdf"]}
    tp_ylog_norm_plotspecs.update(plotstyle)

    # gather data
    tp_data = project_util.Dataset(topologies, top_names, b_values, replicas, datapath)

    # normalize data 
    tp_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            tp_0 = np.array(tp_data.data[t][n][0])
            avg_tp_0 = np.mean(tp_0[~np.isnan(tp_0)])
            for b in range(len(b_values)):
                tp_b = tp_data.data[t][n][b]
                for rep in range(len(tp_b)):
                    if not np.isnan(tp_b[rep]):
                        tp_norm_data.data[t][n][b][rep] = tp_b[rep]/avg_tp_0

    tp_ylog_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            tp_0 = np.array(tp_data.data[t][n][0])
            avg_tp_0 = np.mean(tp_0[~np.isnan(tp_0)])
            for b in range(len(b_values)):
                tp_b = tp_data.data[t][n][b]
                for rep in range(len(tp_b)):
                    if not np.isnan(tp_b[rep]):
                        tp_ylog_norm_data.data[t][n][b][rep] = np.log(tp_b[rep]/avg_tp_0)
    tp_ylog_norm_data._calc_repavg()



    # plot data
    project_plotter.plot_data(tp_ylog_norm_data, tp_ylog_norm_plotspecs)
    project_plotter.plot_data(tp_norm_data, tp_norm_plotspecs)
    project_plotter.plot_data(tp_data, tp_plotspecs)
    plt.show()
