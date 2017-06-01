import os
import numpy as np
import matplotlib.pyplot as plt

import project_plotter
import project_util

if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]
    replicas = range(1,11)

    plotstyle = project_plotter.global_plot_style()

    datapath = "msm/n_physical_timescales.dat"
    n_phys_ti_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"Num. of physical timescales", 
            "title":" ", "legend_loc":2, "ylims":(0,10),
            "saveas":"msm/repavg_msm_n_phys_times", "saveas_formats":["png","pdf"]}
    n_phys_ti_plotspecs.update(plotstyle)
    
    n_phys_ti_data = project_util.Dataset(topologies, top_names, b_values, replicas, datapath)
    project_plotter.plot_data(n_phys_ti_data, n_phys_ti_plotspecs)

    datapath = "msm/n_dominant_timescales.dat"
    n_dom_ti_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"Num. of dominant timescales",
            "title":" ", "legend_loc":2, "ylims":(0,3),
            "saveas":"msm/repavg_msm_n_dom_times", "saveas_formats":["png","pdf"]}

    n_dom_ti_plotspecs.update(plotstyle)

    n_dom_ti_data = project_util.Dataset(topologies, top_names, b_values, replicas, datapath)
    project_plotter.plot_data(n_dom_ti_data, n_dom_ti_plotspecs)

    plt.show()
