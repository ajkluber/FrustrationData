import os
import numpy as np
import matplotlib.pyplot as plt

import project_util, project_plotter

if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]

    replicas = range(1,11)

    plotstyle = project_plotter.global_plot_style()

    msm_mfpt_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$t_{mfpt}$ (frames)", 
            "title":"MFPT along Q", "legend_loc":6,
            "saveas":"repavg_msm_mfpt_1D", "saveas_formats":["png","pdf"]}
    msm_mfpt_plotspecs.update(plotstyle)

    msm_mfpt_ylog_norm_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$\\log\\frac{t_{mfpt}}{t^{b=0}_{mfpt}}$", 
            "title":"MFPT along Q", "legend_loc":2,
            "saveas":"repavg_msm_mfpt_1D_ylog_norm", "saveas_formats":["png","pdf"]}
    msm_mfpt_ylog_norm_plotspecs.update(plotstyle)

    datapath = "msm/tfold_pcca.dat"

    msm_mfpt_data = project_util.Dataset(topologies, top_names, b_values, replicas, datapath)

    msm_mfpt_ylog_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            mfpt_0 = msm_mfpt_data.data[t][n][0]
            avg_mfpt_0 = np.mean(mfpt_0[~np.isnan(mfpt_0)])
            for b in range(len(b_values)):
                mfpt_b = msm_mfpt_data.data[t][n][b]
                for rep in range(len(mfpt_b)):
                    y = mfpt_b[rep]
                    if np.isnan(y):
                        msm_mfpt_ylog_norm_data.data[t][n][b][rep] = np.nan
                    else:
                        msm_mfpt_ylog_norm_data.data[t][n][b][rep] = np.log(y/avg_mfpt_0)
    msm_mfpt_ylog_norm_data._calc_repavg()

    project_plotter.plot_data(msm_mfpt_ylog_norm_data, msm_mfpt_ylog_norm_plotspecs)
    project_plotter.plot_data(msm_mfpt_data, msm_mfpt_plotspecs)
