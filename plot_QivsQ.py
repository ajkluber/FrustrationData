import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid.inset_locator import inset_axes

import project_util
import project_plotter

class DatasetXvsY(project_util.DatasetXvsY):
    def __init__(self, *args, **kwargs):
        project_util.DatasetXvsY.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        data_x = np.loadtxt(self.datapath_x)
        data_y = np.load(self.datapath_y)
        #M = float(len(np.loadtxt("Qi_vs_Qtanh_0_05/nonnative_pairs.dat")))
        #data_y *= M
        return data_x, data_y

def subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, dataset):
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*4.7, n_rows*3.5))
    for i in range(n_rows):
        for j in range(n_cols):
            t = data_specs["topologies"][i][j] 
            n = data_specs["top_names"][i][j] 
            b = data_specs["b_values"][i][j]

            pairs = dataset.pairs[t][n]
            N = dataset.prot_sizes[t][n]
            #print dataset.top_names[t][n]

            #print dataset.ydata[t][n][b]
            ax = axes[i][j]
            #import pdb; pdb.set_trace()
            print t, n, b
            vals = dataset.ydata[t][n][b][0]
            C = np.zeros((N, N))
            for m in range(len(pairs)):
                if m < dataset.prot_n_native[t][n]:
                    C[pairs[m, 1], pairs[m, 0]] = vals[m]
                else:
                    C[pairs[m, 1], pairs[m, 0]] = -vals[m]
            
            # plot native and non-native contacts in different colors
            vmin, vmax = subplotspecs["vminmax"]
            pa = ax.pcolormesh(np.ma.array(C, mask=(C == 0)), cmap="bwr_r", vmin=vmin, vmax=vmax)

            ax.plot(np.arange(0, N + 1), np.arange(0, N + 1), 'k', lw=2)

            if subplotspecs.has_key("zoom"):
                fact = subplotspecs["zoom_factor"]
                ax_in = inset_axes(ax, "40%", "45%", loc=4)
                x_idxs, y_idxs = subplotspecs["zoom"][i]
                Cinset = C[y_idxs[0]:y_idxs[1] + 1, x_idxs[0]:x_idxs[1] + 1]
                X, Y = np.meshgrid(np.arange(x_idxs[0], x_idxs[1] + 1), np.arange(y_idxs[0], y_idxs[1] + 1))
                p_in = ax_in.pcolormesh(X, Y, fact*Cinset, cmap="bwr_r", vmin=vmin, vmax=vmax)
                ax_in.set_xticks([])
                ax_in.set_yticks([])
                if subplotspecs.has_key("zoom_title"):
                    ax_in.set_title(subplotspecs["zoom_title"], fontsize=16)

                #cbar_in = plt.colorbar(p_in, ax_in, ticks=[])
                #cbar_in.set_label("2x scale")


                ax.add_patch(patches.Rectangle((x_idxs[0], y_idxs[0]), 
                    y_idxs[1] - y_idxs[0], 
                    x_idxs[1] - x_idxs[0], fill=False, lw=1.2, ec="k"))


            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.set_xlim(0, N)
            ax.set_ylim(0, N)

            # add labels
            if j > 0:
                ax.set_yticklabels([])
            else:
                ax.set_ylabel("Residue $j$", fontsize=12)

            if i == (n_rows - 1):
                ax.set_xlabel("Residue $i$", fontsize=12)

            if j == 0: 
                ax.annotate(subplotspecs["panel_labels"][i],
                    xy=(0,0), xytext=(-0.3, 0.925), fontsize=26,
                    xycoords="axes fraction", textcoords="axes fraction")

                ax.annotate(subplotspecs["prot_labels"][i],
                    xy=(0,0), xytext=(-0.34, 0.5), fontsize=26,
                    xycoords="axes fraction", textcoords="axes fraction")

            if i == 0:
                ax.annotate(r"$b = {}$".format(dataset.b_values[b]),
                    xy=(0,0), xytext=(0.4, 1.05), fontsize=20,
                    xycoords="axes fraction", textcoords="axes fraction")

    plt.subplots_adjust(wspace=0.075)
    fig.subplots_adjust(left=0.05, right=0.9)
    cbar_ax1 = fig.add_axes([0.91, 0.1, 0.01, 0.8])
    cb = fig.colorbar(pa, cax=cbar_ax1, ticks=subplotspecs["cbar_ticks"])
    cb.ax.set_yticklabels(subplotspecs["cbar_ticklabels"])
    cbar_ax1.annotate("Native",
                    xy=(0,0), xytext=(2, 0.1 + 0.7), fontsize=16,
                    xycoords="axes fraction", textcoords="axes fraction",
                    rotation="vertical")
    cbar_ax1.annotate("Non-Native",
                    xy=(0,0), xytext=(2, 0.1 + 0.25), fontsize=16,
                    xycoords="axes fraction", textcoords="axes fraction",
                    rotation="vertical")

    for format in [".png", ".pdf"]:
        fig.savefig("plots/" + subplotspecs["saveas"] + format)
        
def plot_1r69_1fmk_maps(TSavg_data, TSstd_data, Uavg_data, Ustd_data):
    # plot subplot of contact maps particular values
    n_rows, n_cols = 2, 3

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[0, 0, 0], [0, 0, 0]],
            "b_values":[[3, 6, 10], [3, 6, 10]]}
    subplotspecs = {"vminmax":(-1, 1),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_1$", r"$\beta_1$"],
            "panel_labels":["A", "B"],
            "cbar_ticks":np.linspace(-1, 1, 3), 
            "cbar_ticklabels":["1", "0", "1"],
            "zoom":(((0,10), (28,38)), ((8,17), (43,52))),
            "zoom_factor":2,
            "zoom_title":"2x scale",
            "saveas":"grid_TSavg_1r69_1fmk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, TSavg_data)

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[0, 0, 0], [0, 0, 0]],
            "b_values":[[3, 6, 10], [3, 6, 10]]}

    subplotspecs = {"vminmax":(-.5, 0.5),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_1$", r"$\beta_1$"],
            "panel_labels":["A", "B"],
            "cbar_ticks":np.linspace(-.5, 0.5, 3),
            "cbar_ticklabels":["0.5", "0", "0.5"],
            "zoom":(((0,10), (28,38)), ((8,17), (43,52))),
            "zoom_factor":1,
            "saveas":"grid_TSstd_1r69_1fmk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, TSstd_data)

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[0, 0, 0], [0, 0, 0]],
            "b_values":[[3, 6, 10], [3, 6, 10]]}
    subplotspecs = {"vminmax":(-1, 1),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_1$", r"$\beta_1$"],
            "panel_labels":["A", "B"],
            "cbar_ticks":np.linspace(-1, 1, 3),
            "cbar_ticklabels":["1", "0", "1"],
            "zoom":(((0,10), (28,38)), ((8,17), (43,52))),
            "zoom_factor":2,
            "zoom_title":"2x scale",
            "saveas":"grid_Uavg_1r69_1fmk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, Uavg_data)

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[0, 0, 0], [0, 0, 0]],
            "b_values":[[3, 6, 10], [3, 6, 10]]}
    subplotspecs = {"vminmax":(-.5, 0.5),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_1$", r"$\beta_1$"],
            "panel_labels":["A", "B"],
            "cbar_ticks":np.linspace(-.5, 0.5, 3),
            "cbar_ticklabels":["0.5", "0", "0.5"],
            "zoom":(((0,10), (28,38)), ((8,17), (43,52))),
            "zoom_factor":1,
            "saveas":"grid_Ustd_1r69_1fmk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, Ustd_data)

def plot_1im1_2akk_maps(TSavg_data, TSstd_data, Uavg_data, Ustd_data):
    # plot subplot of contact maps particular values
    n_rows, n_cols = 2, 3

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[1, 1, 1], [1, 1, 1]],
            "b_values":[[7, 12, 14], [7, 12, 14]]}
    subplotspecs = {"vminmax":(-1, 1),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_2$", r"$\beta_2$"],
            "panel_labels":["A", "B"],
            "cbar_ticks":np.linspace(-1, 1, 3), 
            "cbar_ticklabels":["1", "0", "1"],
            "zoom":(((5,21), (66,82)), ((0,13), (40,53))),
            "zoom_factor":2,
            "zoom_title":"2x scale",
            "saveas":"grid_TSavg_1imq_2akk"}
            #"zoom":(((38,48), (60,70)), ((8,17), (43,52))), # helix 3 and helix 4

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, TSavg_data)

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[1, 1, 1], [1, 1, 1]],
            "b_values":[[7, 12, 14], [7, 12, 14]]}
    subplotspecs = {"vminmax":(-.5, 0.5),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_2$", r"$\beta_2$"],
            "panel_labels":["A", "B"],
            "cbar_ticks":np.linspace(-.5, 0.5, 3),
            "cbar_ticklabels":["0.5", "0", "0.5"],
            "zoom":(((5,21), (66,82)), ((0,13), (40,53))),
            "zoom_factor":1,
            "saveas":"grid_TSstd_1imq_2akk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, TSstd_data)

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[1, 1, 1], [1, 1, 1]],
            "b_values":[[7, 12, 14], [7, 12, 14]]}
    subplotspecs = {"vminmax":(-1, 1),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_2$", r"$\beta_2$"],
            "panel_labels":["A", "B"],
            "cbar_ticks":np.linspace(-1, 1, 3),
            "cbar_ticklabels":["1", "0", "1"],
            "zoom":(((5,21), (66,82)), ((0,13), (40,53))),
            "zoom_factor":2,
            "zoom_title":"2x scale",
            "saveas":"grid_Uavg_1imq_2akk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, Uavg_data)

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[1, 1, 1], [1, 1, 1]],
            "b_values":[[7, 12, 14], [7, 12, 14]]}
    subplotspecs = {"vminmax":(-.5, 0.5),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_2$", r"$\beta_2$"],
            "panel_labels":["A", "B"],
            "cbar_ticks":np.linspace(-.5, 0.5, 3),
            "cbar_ticklabels":["0.5", "0", "0.5"],
            "zoom":(((5,21), (66,82)), ((0,13), (40,53))),
            "zoom_factor":1,
            "saveas":"grid_Ustd_1imq_2akk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, Ustd_data)

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

    plotstyle = project_plotter.global_plot_style()

    replicas = range(1,11)

    coordname = "Qtanh_0_05"
    datapath_x = "Qi_vs_" + coordname + "/Q_mid_bin.dat"
    datapath_y = "Qi_vs_" + coordname + "/Qi_vs_Q.npy"
    grid_dims = (4, 4)
    ylims = (0,1)

    #saveas = []
    #for i in range(len(topologies)):
    #    temp = []
    #    for j in range(len(top_names[i])):
    #        temp.append("repavg_QvsQtanh_{}".format(top_names[i][j]))
    #    saveas.append(temp)

    Q_plotspecs = {"xlabel":"Native Contacts", "ylabel":"$Q_i$", 
            "title":"Folding mechanism", "saveas":"qi_vs_Q_replicas",
            "grid_dims":grid_dims, "coloridxs":coloridxs, "ylims":ylims, "xlims":(0,1),
            "xytext":(0.7,0.05)}

    TScont_plotspecs = {"xlabel":"index", "ylabel":"index", 
            "title":"TS structure", "saveas":"q_TS_replicas",
            "grid_dims":grid_dims, "coloridxs":coloridxs,
            "xytext":(0.6,0.12), "cmap":"Blues", "vminmax":((0, 1), (0, 1))}
    TScont_plotspecs.update(plotstyle)

    TSavg_plotspecs = {"xlabel":"index", "ylabel":"index", 
            "title":"Avg TS", "saveas":"repavg_TScont_{}",
            "grid_dims":grid_dims, "coloridxs":coloridxs,
            "xytext":(0.6,0.12), "cmap":"Blues", "vminmax":(-1, 1)}
    TSavg_plotspecs.update(plotstyle)

    TSstd_plotspecs = {"xlabel":"index", "ylabel":"index", 
            "title":"Std dev TS", "saveas":"repstd_TScont_{}",
            "grid_dims":grid_dims, "coloridxs":coloridxs,
            "xytext":(0.6,0.12), "cmap":"Blues", "vminmax":(0, 0.5)}
    TSstd_plotspecs.update(plotstyle)

    Uavg_plotspecs = {"xlabel":"index", "ylabel":"index", 
            "title":"Avg U", "saveas":"repavg_Ucont_{}",
            "grid_dims":grid_dims, "coloridxs":coloridxs,
            "xytext":(0.6,0.12), "cmap":"Blues", "vminmax":(-1, 1)}
    Uavg_plotspecs.update(plotstyle)

    Ustd_plotspecs = {"xlabel":"index", "ylabel":"index", 
            "title":"Std dev U", "saveas":"repstd_Ucont_{}",
            "grid_dims":grid_dims, "coloridxs":coloridxs,
            "xytext":(0.6,0.12), "cmap":"Blues", "vminmax":(0, 0.5)}
    Ustd_plotspecs.update(plotstyle)

    Qi_raw_data = DatasetXvsY(topologies, top_names, b_values, replicas, datapath_x,
                datapath_y, calc_repavg=False)

    Ai_raw_data = DatasetXvsY(topologies, top_names, b_values, replicas, 
            "Ai_vs_{}/Q_mid_bin.dat".format(coordname), "Ai_vs_{}/Ai_vs_Q.npy".format(coordname),
            calc_repavg=False)

    TScont_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    all_pairs = []
    for t in range(len(topologies)):
        top_pairs = []
        for n in range(len(top_names[t])):
            #print top_names[t][n]
            #print "b-value " + "".join([ "{:8}".format(x) for x in replicas ])
            for b in range(len(b_values)):
                repstring = "{:8}".format(b_values[b]) 
                if b == 0:
                    pairs = np.loadtxt("{0}/{1}/b_{2}/replica_1/{1}_pairwise_params".format(topologies[t], top_names[t][n], b_values[b]), usecols=(0,1), dtype=int) - 1
                    top_pairs.append(pairs)
                for rep in range(len(Qi_raw_data.ydata[t][n][b])):
                    M = float(prot_n_native[t][n])

                    # for simplicity take 'transition state' as the Q value
                    # half-way between the outermost minima.
                    temp = np.loadtxt("{}/{}/b_{}/replica_{}/{}_profile/minima.dat".format(
                        topologies[t], top_names[t][n], b_values[b], rep + 1, coordname))
                    left_min, right_min = temp.min(), temp.max()
                    middle = 0.5*(left_min + right_min)
                    repstring += "{:>8.4f}".format(middle/M)

                    qTS_idx = np.argmin((Qi_raw_data.xdata[t][n][b][rep] - middle)**2)
                    aTS_idx = np.argmin((Ai_raw_data.xdata[t][n][b][rep] - middle)**2)
                    qi_TS = Qi_raw_data.ydata[t][n][b][rep][qTS_idx]
                    Ai_TS = Ai_raw_data.ydata[t][n][b][rep][aTS_idx]
                    TScont_data.ydata[t][n][b].append(np.concatenate((qi_TS, Ai_TS)))
                #print repstring
        all_pairs.append(top_pairs)
    TScont_data.pairs = all_pairs
    TScont_data.prot_sizes = prot_sizes
    TScont_data.prot_n_native = prot_n_native

    # calculate average TS contact map and std dev.
    TSavg_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    TSstd_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    all_pairs = []
    for t in range(len(topologies)):
        top_pairs = []
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                if b == 0:
                    pairs = np.loadtxt("{0}/{1}/b_{2}/replica_1/{1}_pairwise_params".format(topologies[t], top_names[t][n], b_values[b]), usecols=(0,1), dtype=int) - 1
                    top_pairs.append(pairs)

                # calculate average mechanism and fluctuations around mechanism
                if len(TScont_data.ydata[t][n][b]) > 0:
                    vals = np.array(TScont_data.ydata[t][n][b])
                    TSavg_data.ydata[t][n][b].append(np.mean(vals, axis=0)) 
                    TSstd_data.ydata[t][n][b].append(np.std(vals, axis=0)) 

        all_pairs.append(top_pairs)
    TSavg_data.pairs = all_pairs
    TSavg_data.prot_sizes = prot_sizes
    TSavg_data.prot_n_native = prot_n_native
    TSstd_data.pairs = all_pairs
    TSstd_data.prot_sizes = prot_sizes
    TSstd_data.prot_n_native = prot_n_native

    # calculate average and fluctuations for unfolded state.
    Ucont_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    all_pairs = []
    for t in range(len(topologies)):
        top_pairs = []
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                if b == 0:
                    pairs = np.loadtxt("{0}/{1}/b_{2}/replica_1/{1}_pairwise_params".format(topologies[t], top_names[t][n], b_values[b]), usecols=(0,1), dtype=int) - 1
                    top_pairs.append(pairs)
                for rep in range(len(Qi_raw_data.ydata[t][n][b])):
                    M = float(prot_n_native[t][n])

                    # for simplicity take 'transition state' as the Q value
                    # half-way between the outermost minima.
                    temp = np.loadtxt("{}/{}/b_{}/replica_{}/{}_profile/minima.dat".format(
                        topologies[t], top_names[t][n], b_values[b], rep + 1, coordname))
                    left = temp.min()

                    qU_idx = np.argmin((Qi_raw_data.xdata[t][n][b][rep] - left)**2)
                    aU_idx = np.argmin((Ai_raw_data.xdata[t][n][b][rep] - left)**2)
                    qi_U = Qi_raw_data.ydata[t][n][b][rep][qU_idx]
                    Ai_U = Ai_raw_data.ydata[t][n][b][rep][aU_idx]
                    Ucont_data.ydata[t][n][b].append(np.concatenate((qi_U, Ai_U)))
        all_pairs.append(top_pairs)
    Ucont_data.pairs = all_pairs
    Ucont_data.prot_sizes = prot_sizes
    Ucont_data.prot_n_native = prot_n_native

    # calculate average U contact map and std dev.
    Uavg_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    Ustd_data = project_util.DatasetXvsY(topologies, top_names, b_values, replicas,"","", empty=True)
    all_pairs = []
    for t in range(len(topologies)):
        top_pairs = []
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                if b == 0:
                    pairs = np.loadtxt("{0}/{1}/b_{2}/replica_1/{1}_pairwise_params".format(topologies[t], top_names[t][n], b_values[b]), usecols=(0,1), dtype=int) - 1
                    top_pairs.append(pairs)

                # calculate average mechanism and fluctuations around mechanism
                if len(Ucont_data.ydata[t][n][b]) > 0:
                    vals = np.array(Ucont_data.ydata[t][n][b])
                    Uavg_data.ydata[t][n][b].append(np.mean(vals, axis=0)) 
                    Ustd_data.ydata[t][n][b].append(np.std(vals, axis=0)) 

        all_pairs.append(top_pairs)
    Uavg_data.pairs = all_pairs
    Uavg_data.prot_sizes = prot_sizes
    Uavg_data.prot_n_native = prot_n_native
    Ustd_data.pairs = all_pairs
    Ustd_data.prot_sizes = prot_sizes
    Ustd_data.prot_n_native = prot_n_native

    # plot fancy figures for paper
    #plot_1im1_2akk_maps(TSavg_data, TSstd_data, Uavg_data, Ustd_data)
    #plot_1r69_1fmk_maps(TSavg_data, TSstd_data, Uavg_data, Ustd_data)

    # calculate route measure

    #project_plotter.plot_bvalue_maps_grid(TScont_data, Q_TS_plotspecs)

    project_plotter.plot_replica_maps_grid(TScont_data, TScont_plotspecs) # plot contact maps for all replicas

    #project_plotter.plot_bvalue_maps_grid(TSavg_data, TSavg_plotspecs) # plot avg contact map
    #project_plotter.plot_bvalue_maps_grid(TSstd_data, TSstd_plotspecs) # plot avg contact map
    #project_plotter.plot_bvalue_maps_grid(Uavg_data, Uavg_plotspecs) # plot avg contact map
    #project_plotter.plot_bvalue_maps_grid(Ustd_data, Ustd_plotspecs) # plot std contact map
    #plt.show()

