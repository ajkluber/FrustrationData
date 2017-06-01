import os
import glob
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter

class Dataset(project_util.Dataset):
    def __init__(self, *args, **kwargs):
        project_util.Dataset.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        with open(self.datapath, "r") as fin:
            data = abs(float(fin.read()))
        return data

if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    markers = ["^", "o", "s"]
    N = np.array([63., 86., 58., 74., 48.])

    plotstyle = project_plotter.global_plot_style()

    replicas = range(1,11)
    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    corr_data = Dataset(topologies, top_names, b_values, replicas, "msm/q_vs_psi1_corr.dat")
    corr_TS_data = Dataset(topologies, top_names, b_values, replicas, "msm/q_vs_psi1_corr_TS.dat")

    # Plot correlation of Q with 
    corr_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$\langle Q \cdot\psi_2\rangle$", "ylims": (0.5,1),
            "title":" ", "legend_loc":3,
            "saveas":"corr_Q_vs_psi2", "saveas_formats":["png","pdf"]}
    corr_plotspecs.update(plotstyle)

    corr_TS_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$\langle Q \cdot\psi_2\rangle_{TS}$", "ylims": (0.5,1),
            "title":" ", "legend_loc":3,
            "saveas":"corr_Q_vs_psi2_TS", "saveas_formats":["png","pdf"]}
    corr_TS_plotspecs.update(plotstyle)

    #project_plotter.plot_data(corr_data, corr_plotspecs)
    #project_plotter.plot_data(corr_TS_data, corr_TS_plotspecs)
    #plt.show()

    corr_TS_plotspecs.pop("legend_key")
    corr_TS_plotspecs.pop("legend_loc")

    both_datasets = [corr_data, corr_TS_data]
    both_plotspecs = [corr_plotspecs, corr_TS_plotspecs]

    panel_label = ["A", "B"]

    # plot two colm panel
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for d in [0,1]: 
        dataset = both_datasets[d]
        plotspecs = both_plotspecs[d]
        ax = axes[d]
        for t in range(len(dataset.topologies)):
            names = dataset.top_names[t]
            for n in range(len(names)):
                first = True
                b_line = []
                avg_val_line = []
                std_dev_line = []
                for j in range(len(dataset.b_values)):
                    b = float(dataset.b_values[j])

                    if hasattr(dataset, "pre_averaged"):
                        avg_val = dataset.data[t][n][j]
                        std_val = 0
                    else:
                        if len(dataset.data[t][n][j]) > 0:
                            rep_data = dataset.data[t][n][j]
                            avg_val = np.mean(rep_data)
                            std_val = np.std(rep_data)/np.sqrt(len(rep_data))
                        else:
                            print "no values for b = ", b
                            avg_val = None

                    # plot average value with error bars 
                    if not (avg_val is None):
                        b_line.append(b)
                        avg_val_line.append(avg_val)
                        std_dev_line.append(std_val)

                if plotspecs.has_key("legend_key"):
                    ax.errorbar(b_line, avg_val_line, yerr=std_dev_line,
                            color=plotspecs["color"][t][n],
                            marker=plotspecs["markers"][t], 
                            label=plotspecs["legend_key"][t][n], lw=2)
                else:
                    ax.errorbar(b_line, avg_val_line, yerr=std_dev_line,
                            color=plotspecs["color"][t][n],
                            marker=plotspecs["markers"][t], lw=2)

        #ax.set_ylim(-1.5, 1.5)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 

        if plotspecs.has_key("legend_loc"):
            ax.legend(loc=plotspecs["legend_loc"], fontsize=16)

        ax.set_ylim(plotspecs["ylims"])
        ax.set_ylabel(plotspecs["ylabel"], fontsize=16)
        ax.set_xlabel(plotspecs["xlabel"], fontsize=16)
        plt.subplots_adjust(right=0.98, left=0.02, wspace=0.25)
        
        ax.annotate(panel_label[d], xy=(0,0), xytext=(-0.17, 1.05), fontsize=26,
            xycoords="axes fraction", textcoords="axes fraction")
    #plt.tight_layout()

    for format in plotspecs["saveas_formats"]:
        fig.savefig("plots/two_colm_corrQvsPsi2" + "." + format, bbox_inches="tight")
    plt.show()
    
