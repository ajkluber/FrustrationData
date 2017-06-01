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
        return data

    def _calc_repavg(self):
        """ calculate the replicas average data"""

        all_data_avg = []
        all_data_std = []
        for t in range(len(self.topologies)):
            data_top_avg = []
            data_top_std = []
            for n in range(len(self.top_names[t])):
                data_name_avg = []
                data_name_std = []
                for b in range(len(self.b_values)):
                    # take longest
                    if len(self.data[t][n][b]) > 0:
                        n_reps = len(self.data[t][n][b])
                        maxL = int(np.max([ len(x) for x in self.data[t][n][b] ]))
                        fulldata = np.ones((n_reps, maxL))
                        for rep in range(len(self.data[t][n][b])):
                            if len(self.data[t][n][b][rep]) == 1:
                                pass
                            else:
                                fulldata[rep,:len(self.data[t][n][b][rep])] = self.data[t][n][b][rep]

                        data_name_avg.append(np.mean(fulldata, axis=0))
                        data_name_std.append(np.std(fulldata, axis=0)/np.sqrt(float(n_reps)))
                    else:
                        data_name_avg.append([])
                        data_name_std.append([])

                data_top_avg.append(data_name_avg)
                data_top_std.append(data_name_std)

            all_data_avg.append(data_top_avg)
            all_data_std.append(data_top_std)

        self.avgdata = all_data_avg
        self.stddata = all_data_std

if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]
    replicas = range(1,11)

    plotstyle = project_plotter.global_plot_style()

    datapath = "msm/KC.npy"
    KC_plotspecs = {"xlabel":"Index", "ylabel":"KC",
            "legend_loc":4, "ylims":(0,1), "saveas":"msm/repavg_KC"}
    KC_plotspecs.update(plotstyle)
    
    KC_data = Dataset(topologies, top_names, b_values, replicas, datapath)

    fig, axes = plt.subplots(2,3, figsize=(18,6))
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            ax = axes[n][t]
            for b in range(len(b_values)):
                if len(KC_data.avgdata[t][n][b]) > 0:
                    if b == 0:
                        ax.plot(KC_data.avgdata[t][n][b], color=plotstyle["color"][t][n], label=KC_plotspecs["legend_key"][t][n])
                    ax.plot(KC_data.avgdata[t][n][b], color=plotstyle["color"][t][n])
                    x = np.arange(len(KC_data.avgdata[t][n][b]))
                    y1 = KC_data.avgdata[t][n][b] - KC_data.stddata[t][n][b]
                    y2 = KC_data.avgdata[t][n][b] + KC_data.stddata[t][n][b]
                    ax.fill_between(x, y1, y2, color=plotstyle["color"][t][n], alpha=0.2)
                    ax.set_xlim(0, 25)
                    ax.set_ylim(0, 1)
            ax.annotate(KC_plotspecs["legend_key"][t][n],
                    xy=(0,0), xytext=(0.5, 0.1), fontsize=24,
                    xycoords="axes fraction", textcoords="axes fraction")

            if t == 0:
                ax.set_ylabel(KC_plotspecs["ylabel"], fontsize=20)
            if n == 1:
                ax.set_xlabel(KC_plotspecs["xlabel"], fontsize=20)
    fig.suptitle("Kinetic content", fontsize=20)
    fig.savefig("plots/KC_grid.pdf")
    fig.savefig("plots/KC_grid.png")
    plt.show()

    #ls = ["-", "--"]
    #for t,n in [(0,0), (1,0)]:
    #    for b in [3,6,10]:
    #    for b in range(len(b_values)):
    #        if b == 0:
    #            plt.plot(KC_data.avgdata[t][n][b], ls=ls[t], color=plotstyle["color"][t][n], label=KC_plotspecs["legend_key"][t][n])
    #        else:
    #            plt.plot(KC_data.avgdata[t][n][b], ls=ls[t], color=plotstyle["color"][t][n])
    #        plt.xlim(0, 25)
    #        plt.ylim(0, 1)
    #plt.legend(loc=4)
    #plt.xlabel(KC_plotspecs["xlabel"])
    #plt.ylabel(KC_plotspecs["ylabel"])
    #plt.show()
    #project_plotter.plot_data(KC_data, KC_plotspecs)
    #plt.show()
