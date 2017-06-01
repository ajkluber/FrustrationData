import os
import numpy as np
import matplotlib.pyplot as plt

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

            dataset.ydata[t][n][b]
            ax = axes[i][j]

            vals = dataset.ydata[t][n][b][0]
            Cnat = np.zeros((N, N))
            Cnon = np.zeros((N, N))
            for m in range(len(pairs)):
                if m < dataset.prot_n_native[t][n]:
                    Cnat[pairs[m, 1], pairs[m, 0]] = vals[m]
                else:
                    Cnon[pairs[m, 1], pairs[m, 0]] = vals[m]
            
            # plot native and non-native contacts in different colors
            vmin, vmax = subplotspecs["vminmax"][0]
            pa = ax.pcolormesh(np.ma.array(Cnat, mask=Cnat == 0), cmap="Blues", vmin=vmin, vmax=vmax)

            vmin, vmax = subplotspecs["vminmax"][1]
            pb = ax.pcolormesh(np.ma.array(Cnon, mask=Cnon == 0), cmap="Reds", vmin=vmin, vmax=vmax)

            ax.plot(np.arange(0, N + 1), np.arange(0, N + 1), 'k', lw=2)
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
    fig.subplots_adjust(left=0.05, right=0.8)
    cbar_ax1 = fig.add_axes([0.82, 0.1, 0.02, 0.8])
    cbar_ax2 = fig.add_axes([0.92, 0.1, 0.02, 0.8])
    cba = fig.colorbar(pa, cax=cbar_ax1, ticks=subplotspecs["cbar1_ticks"])
    cbb = fig.colorbar(pb, cax=cbar_ax2, ticks=subplotspecs["cbar2_ticks"])
    cba.set_label("Native", fontsize=16)
    cbb.set_label("Non-native", fontsize=16)
    for format in [".png", ".pdf"]:
        fig.savefig("plots/" + subplotspecs["saveas"] + format)
        

if __name__ == "__main__":
    topologies = ["alpha","beta", "mixed"]
    top_names = [["1r69", "1imq"],["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    prot_n_native = [[155, 219], [164, 213], [117]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]
    coloridxs = [ b/1.3 for b in float_b ]

    #plotstyle = project_plotter.global_plot_style()

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
            "xytext":(0.6,0.12), "cmap":"Blues", "vminmax":((0, 1), (0, 0.25))}

    TSavg_plotspecs = {"xlabel":"index", "ylabel":"index", 
            "title":"Avg TS", "saveas":"repavg_TScont_{}",
            "grid_dims":grid_dims, "coloridxs":coloridxs,
            "xytext":(0.6,0.12), "cmap":"Blues", "vminmax":((0, 1), (0, 0.25))}

    TSstd_plotspecs = {"xlabel":"index", "ylabel":"index", 
            "title":"Std dev TS", "saveas":"repstd_TScont_{}",
            "grid_dims":grid_dims, "coloridxs":coloridxs,
            "xytext":(0.6,0.12), "cmap":"Blues", "vminmax":((0, 0.5), (0, 0.5))}

    Uavg_plotspecs = {"xlabel":"index", "ylabel":"index", 
            "title":"Avg U", "saveas":"repavg_Ucont_{}",
            "grid_dims":grid_dims, "coloridxs":coloridxs,
            "xytext":(0.6,0.12), "cmap":"Blues", "vminmax":((0, 1), (0, 0.25))}

    Ustd_plotspecs = {"xlabel":"index", "ylabel":"index", 
            "title":"Std dev U", "saveas":"repstd_Ucont_{}",
            "grid_dims":grid_dims, "coloridxs":coloridxs,
            "xytext":(0.6,0.12), "cmap":"Blues", "vminmax":((0, 0.5), (0, 0.5))}

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
            print top_names[t][n]
            print "b-value " + "".join([ "{:8}".format(x) for x in replicas ])
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
                print repstring
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

    # plot subplot of contact maps particular values
    n_rows, n_cols = 2, 3

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[0, 0, 0], [0, 0, 0]],
            "b_values":[[3, 6, 10], [3, 6, 10]]}
    subplotspecs = {"vminmax":((0, 1), (0, 0.25)),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_1$", r"$\beta_1$"],
            "panel_labels":["A", "B"],
            "cbar1_ticks":np.arange(0, 1.2, 0.2),
            "cbar2_ticks":np.arange(0, 0.25 + 0.05, 0.05),
            "saveas":"grid_TSavg_1r69_1fmk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, TSavg_data)

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[0, 0, 0], [0, 0, 0]],
            "b_values":[[3, 6, 10], [3, 6, 10]]}

    subplotspecs = {"vminmax":((0, 0.5), (0, 0.5)),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_1$", r"$\beta_1$"],
            "panel_labels":["A", "B"],
            "cbar1_ticks":np.arange(0, 0.6, 0.1),
            "cbar2_ticks":np.arange(0, 0.6, 0.1),
            "saveas":"grid_TSstd_1r69_1fmk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, TSstd_data)

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[0, 0, 0], [0, 0, 0]],
            "b_values":[[3, 6, 10], [3, 6, 10]]}
    subplotspecs = {"vminmax":((0, 1), (0, 0.25)),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_1$", r"$\beta_1$"],
            "panel_labels":["A", "B"],
            "cbar1_ticks":np.arange(0, 1.2, 0.2),
            "cbar2_ticks":np.arange(0, 0.25 + 0.05, 0.05),
            "saveas":"grid_Uavg_1r69_1fmk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, Uavg_data)

    data_specs = {"topologies":[[0, 0, 0], [1, 1, 1]],
            "top_names":[[0, 0, 0], [0, 0, 0]],
            "b_values":[[3, 6, 10], [3, 6, 10]]}
    subplotspecs = {"vminmax":((0, 1), (0, 0.25)),
            "xytext":(0.7,0.05), "prot_labels":[r"$\alpha_1$", r"$\beta_1$"],
            "panel_labels":["A", "B"],
            "cbar1_ticks":np.arange(0, 1.2, 0.2),
            "cbar2_ticks":np.arange(0, 0.25 + 0.05, 0.05),
            "saveas":"grid_Ustd_1r69_1fmk"}

    subplot_contact_maps(n_rows, n_cols, data_specs, subplotspecs, Ustd_data)
    plt.show()

    # calculate route measure

    #project_plotter.plot_bvalue_maps_grid(TScont_data, Q_TS_plotspecs)

    #plot_replica_maps_grid(TScont_data, TScont_plotspecs) # plot contact maps for all replicas

    #plot_bvalue_maps_grid(TSavg_data, TSavg_plotspecs) # plot avg contact map
    #plot_bvalue_maps_grid(TSavg_data, TSavg_plotspecs) # plot avg contact map
    #plot_bvalue_maps_grid(Uavg_data, Uavg_plotspecs) # plot avg contact map
    #plot_bvalue_maps_grid(Ustd_data, Ustd_plotspecs) # plot std contact map
    #plt.show()

