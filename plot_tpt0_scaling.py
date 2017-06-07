import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

import project_util
import project_plotter

if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    N = np.array([63., 86., 58., 74., 48.])

    plotstyle = project_plotter.global_plot_style()

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25", 
                "1.35", "1.45"]

    replicas = range(1,11)

    tp_data = project_util.Dataset(topologies, top_names, b_values, 
            replicas, "Qtanh_0_05_transit_time/forward_mean")

    # fit scaling law
    avg_unfrust_tp = []
    std_unfrust_tp = []
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            avg_unfrust_tp.append(np.mean(tp_data.data[t][n][0]))
            std_unfrust_tp.append(np.std(tp_data.data[t][n][0])/np.sqrt(len(tp_data.data[t][n][0])))
    avg_unfrust_tp = np.array(avg_unfrust_tp)
    std_unfrust_tp = np.array(std_unfrust_tp)

    line = lambda x, a, b: a + b*x
    popt, pcov = scipy.optimize.curve_fit(line, np.log(N), np.log(avg_unfrust_tp), p0=(-1,0.5))
    a, b = popt

    idx = 0
    plt.figure()
    x = np.linspace(3.8, 4.6, 10)
    plt.plot(x, a + b*x, 'k--')
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            plt.plot(np.log(N[idx]), np.log(avg_unfrust_tp[idx]), 
                    marker=plotstyle["markers"][t], ms=10, color=plotstyle["color"][t][n], 
                    label=plotstyle["legend_key"][t][n])
            idx += 1

    plt.xlabel(r"$\log (N)$")
    plt.ylabel(r"$\log (t_{tp})$")
    plt.title(r"Scaling law  $t_{{tp}} = ({:.2f})$ $N^{{{:.2f}}}$".format(np.exp(a),b))
    plt.legend(loc=2)
    plt.savefig("plots/tp_unfrust_scaling.pdf")
    plt.savefig("plots/tp_unfrust_scaling.png")

    plt.show()
