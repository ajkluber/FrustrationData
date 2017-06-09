import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from plotter.cube_cmap import cubecmap

import matplotlib.cm as cm
cmap = cm.get_cmap("viridis")

def global_plot_style():
    """ """
    # To Do: This function can be smarter and create styles according
    # to a subset of proteins if desired.
    topologies = ["alpha", "beta", "mixed"]
    markers = ["^", "o", "s"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    labels = [[r"$\alpha_1$", r"$\alpha_2$"], [r"$\beta_1$", r"$\beta_2$"], [r"$\alpha/\beta$"]]
    color_codes = [["#e41a1c", "#377eb8"], ["#4daf4a", "#984ea3"], ["#ff7f00"]]

    legend_key = []
    for t in range(len(topologies)):
        top_keys = []
        for n in range(len(top_names[t])):
           #top_keys.append("{} ({})".format(top_names[t][n], prot_sizes[t][n]))
           top_keys.append(labels[t][n])
        legend_key.append(top_keys)

    plotstyle = {"linestyle":"-", "legend_key":legend_key, 
            "markers":markers, "color":color_codes, 
            "newfig":True, "saveas_formats":["png","pdf"]}
    return plotstyle

def plot_repavg(dataset, plotspecs):
    """Plot replica average of obseravables as a function of frustration in a
    grid"""

    grid_dims = plotspecs["grid_dims"]
    for t in range(len(dataset.topologies)):
        names = dataset.top_names[t]
        for n in range(len(names)):
            # Plot whatever for a protein
            fig, axes = plt.subplots(*grid_dims, sharex=True, sharey=True, figsize=(12,10))
            for j in range(len(dataset.b_values)):
                ax = axes[j / grid_dims[1], j % grid_dims[1]]
                rep_xdata = dataset.xdata[t][n][j]
                rep_ydata = dataset.ydata[t][n][j]
                if len(rep_xdata) > 0:
                    # 
                    if hasattr(rep_ydata[0], "mask"):
                        # plot profile for each sample 
                        for r in range(len(rep_xdata)): 
                            #ax.plot(rep_xdata[r], rep_ydata[r], color=cubecmap(plotspecs["coloridxs"][j]), alpha=0.8)
                            ax.plot(rep_xdata[r][~ rep_ydata[r].mask], rep_ydata[r][~ rep_ydata[r].mask], color=cmap(plotspecs["coloridxs"][j]), alpha=0.8)
                    else:
                        # plot profile for each sample 
                        for r in range(len(rep_xdata)): 
                            #ax.plot(rep_xdata[r], rep_ydata[r], color=cubecmap(plotspecs["coloridxs"][j]), alpha=0.8)
                            ax.plot(rep_xdata[r], rep_ydata[r], color=cmap(plotspecs["coloridxs"][j]), alpha=0.8)

                    x_repavg = dataset.avgxdata[t][n][j]
                    y_repavg = dataset.avgydata[t][n][j]
                    # plot disorder-averaged profile
                    ax.plot(x_repavg, y_repavg, lw=2, color='k')

                    if plotspecs.has_key("xytext"):
                        xytext = plotspecs["xytext"]
                    else:
                        xytext = (0.3, 0.1)

                    ax.annotate("b = " + dataset.b_values[j], xy=(0,0), xytext=xytext,
                        bbox={"boxstyle":"square","facecolor":"w","edgecolor":"k"},
                        xycoords="axes fraction", textcoords="axes fraction")

                    if plotspecs.has_key("ylims"): 
                        ax.set_ylim(*plotspecs["ylims"])
                    if plotspecs.has_key("xlims"): 
                        ax.set_xlim(*plotspecs["xlims"])
                    if plotspecs.has_key("ylog"):
                        if plotspecs["ylog"]:
                            ax.semilogy()

            big_ax = fig.add_subplot(111)
            big_ax.grid(False)
            big_ax.set_axis_bgcolor('none')
            big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
            big_ax.set_ylabel(plotspecs["ylabel"])
            big_ax.set_xlabel(plotspecs["xlabel"])
            if plotspecs.has_key("title"):
                big_ax.set_title(plotspecs["title"].format(names[n]))

            plt.subplots_adjust(wspace=0, hspace=0)

            if not (plotspecs["saveas"] is None):
                if not os.path.exists("plots"):
                    os.mkdir("plots")
                os.chdir("plots")
                for format in plotspecs["saveas_formats"]:
                    fig.savefig(plotspecs["saveas"][t][n] + "." + format, bbox_inches="tight")
                os.chdir("..")

             
            #plt.figure()
            #levels = np.linspace(0, 1.3, 11)
            #CS3 = plt.contourf([[0,0],[0,0]], levels, cmap='viridis')
            #plt.clf()

            #plt.figure()
            #for j in range(len(name_repavg_x)):
            #    plt.plot(name_repavg_x[j], name_repavg_y[j], color=cmap(plotspecs["coloridxs"][j]))
            #plt.xlabel(plotspecs["xlabel"])
            #plt.ylabel(plotspecs["ylabel"])
            #plt.title(plotspecs["title"].format(names[n]))

            #cbar = plt.colorbar(CS3)
            #cbar.set_label("Frustration b")
            if plotspecs.has_key("avg_ylims"): 
                plt.ylim(*plotspecs["avg_ylims"])
            if plotspecs.has_key("avg_xlims"): 
                plt.xlim(*plotspecs["avg_xlims"])

            if not (plotspecs["saveas"] is None):
                if not os.path.exists("plots"):
                    os.mkdir("plots")
                os.chdir("plots")
                for format in plotspecs["saveas_formats"]:
                    plt.savefig(plotspecs["saveas"][t][n] + "_avg." + format, bbox_inches="tight")
                os.chdir("..")


def plot_repavg_protein_grid(dataset, plotspecs):
    grid_dims = plotspecs["grid_dims"]
    fig, axes = plt.subplots(*grid_dims, sharex=True, sharey=True, figsize=(12,10))
    prot_idx = 0
    for t in range(len(dataset.topologies)):
        names = dataset.top_names[t]
        for n in range(len(names)):
            # Plot whatever for a protein
            ax = axes[prot_idx % grid_dims[0], prot_idx / grid_dims[0]]
            for j in range(len(dataset.b_values)):
                x_repavg = dataset.avgxdata[t][n][j]
                y_repavg = dataset.avgydata[t][n][j]
                if len(x_repavg) > 0:
                    # plot disorder-averaged profile
                    ax.plot(x_repavg, y_repavg, lw=2, color=cmap(plotspecs["coloridxs"][j]))
            prot_idx += 1

            ax.annotate(names[n] + " (" + dataset.topologies[t] + ")",
                xy=(0,0), xytext=plotspecs["xytext"],
                bbox={"boxstyle":"square","facecolor":"w","edgecolor":"k"},
                xycoords="axes fraction", textcoords="axes fraction")

    plt.subplots_adjust(wspace=0, hspace=0)

    if plotspecs.has_key("ylims"): 
        plt.ylim(*plotspecs["ylims"])
    if plotspecs.has_key("xlims"): 
        plt.xlim(*plotspecs["xlims"])

    big_ax = fig.add_subplot(111)
    big_ax.grid(False)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    big_ax.set_ylabel(plotspecs["ylabel"])
    big_ax.set_xlabel(plotspecs["xlabel"])
    big_ax.set_title(plotspecs["title"])


    if not (plotspecs["saveas"] is None):
        if not os.path.exists("plots"):
            os.mkdir("plots")
        os.chdir("plots")
        for format in plotspecs["saveas_formats"]:
            plt.savefig(plotspecs["saveas"] + "." + format, bbox_inches="tight")
        os.chdir("..")


def plot_data(dataset, plotspecs, err_mult=1.):
    """Plot values versus frustration"""

    if plotspecs["newfig"]:
        fig = plt.figure()
    else:
        fig = plt.gcf()

    for t in range(len(dataset.topologies)):
        names = dataset.top_names[t]
        for n in range(len(names)):
            first = True
            b_line = []
            avg_val_line = []
            std_dev_line = []
            for j in range(len(dataset.b_values)):
                b = float(dataset.b_values[j])

                if hasattr(dataset, "has_avg"):
                    avg_val = dataset.avgdata[t][n][j]
                    std_val = dataset.stddata[t][n][j]
                else:
                    y = dataset.data[t][n][j]
                    good = ~np.isnan(y)
                    if np.sum(good) > 0:
                        avg_val = np.mean(y[good])
                        std_val = np.std(y[good])/np.sqrt(float(np.sum(good)))
                    else:
                        avg_val = np.nan

                # plot average value with error bars 
                if not np.isnan(avg_val):
                    b_line.append(b)
                    avg_val_line.append(avg_val)
                    std_dev_line.append(std_val)
                    if first and plotspecs.has_key("legend_key"):
                        plt.errorbar(b, avg_val, yerr=std_val*err_mult,
                                color=plotspecs["color"][t][n],
                                marker=plotspecs["markers"][t],
                                markersize=8, 
                                label=plotspecs["legend_key"][t][n])
                        first = False
                    else:
                        pass
            
            b_line = np.array(b_line)
            avg_val_line = np.array(avg_val_line)
            std_dev_line = np.array(std_dev_line)
            plt.errorbar(b_line, avg_val_line, yerr=std_dev_line*err_mult,
                        color=plotspecs["color"][t][n],
                        marker=plotspecs["markers"][t], 
                        markersize=8, ls=plotspecs["linestyle"],
                        lw=2)

    if plotspecs.has_key("ylims"): 
        plt.ylim(*plotspecs["ylims"])
    if plotspecs.has_key("xlims"): 
        plt.xlim(*plotspecs["xlims"])

    if plotspecs.has_key("ylog"):
        if plotspecs["ylog"]:
            plt.semilogy()

    if plotspecs.has_key("legend_key"):
        if plotspecs.has_key("legend_loc"):
            plt.legend(loc=plotspecs["legend_loc"], fontsize=18)
        else:
            plt.legend(fontsize=18)

    #plt.text(-0.26, 0.68, plotspecs["ylabel"], fontsize=28) # use for y-label of Rg_packing plot
    #plt.text(-0.26, 0.53, plotspecs["ylabel"], fontsize=30) # use for y-label of packing plot
    plt.ylabel(plotspecs["ylabel"])
    plt.xlabel(plotspecs["xlabel"])
    if plotspecs.has_key("title"):
        plt.title(plotspecs["title"])

    if not (plotspecs["saveas"] is None):
        if not os.path.exists("plots"):
            os.mkdir("plots")
        os.chdir("plots")
        for format in plotspecs["saveas_formats"]:
            fig.savefig(plotspecs["saveas"] + "." + format, bbox_inches="tight")
        os.chdir("..")

def plot_bvalue_grid(dataset, plotspecs):
    """Plot line plots for all replicas of each b-value"""
    cwd = os.getcwd()
    grid_dims = plotspecs["grid_dims"]
    for t in range(len(dataset.topologies)):
        names = dataset.top_names[t]
        for n in range(len(names)):
            # Plot whatever for a protein
            for j in range(len(dataset.b_values)):
                fig, axes = plt.subplots(*grid_dims, sharex=True, sharey=True, figsize=(12,10))
                for rep in range(len(dataset.xdata[t][n][j])):
                    ax = axes[rep / grid_dims[0], rep % grid_dims[0]]
                    # plot
                    x_vals = dataset.xdata[t][n][j][rep]
                    y_vals = dataset.ydata[t][n][j][rep]
                    if len(x_vals) > 0:
                        # plot data for replica
                        ax.plot(x_vals, y_vals)

                    ax.annotate("rep: " + str(rep + 1),
                        xy=(0,0), xytext=plotspecs["xytext"],
                        bbox={"boxstyle":"square","facecolor":"w","edgecolor":"k"},
                        xycoords="axes fraction", textcoords="axes fraction")
                ax.set_aspect("equal")

                if plotspecs.has_key("ylims"): 
                    plt.ylim(*plotspecs["ylims"])
                if plotspecs.has_key("xlims"): 
                    plt.xlim(*plotspecs["xlims"])

                plt.subplots_adjust(wspace=0, hspace=0)
                big_ax = fig.add_subplot(111)
                big_ax.grid(False)
                big_ax.set_axis_bgcolor('none')
                big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
                big_ax.set_ylabel(plotspecs["ylabel"])
                big_ax.set_xlabel(plotspecs["xlabel"])
                big_ax.set_title(plotspecs["title"])

                if not (plotspecs["saveas"] is None):
                    savedir = "{}/{}/b_{}/plots".format(dataset.topologies[t], 
                            dataset.top_names[t][n], dataset.b_values[j])
                    
                    if not os.path.exists(savedir):
                        os.mkdir(savedir)
                    os.chdir(savedir)
                    for format in plotspecs["saveas_formats"]:
                        plt.savefig(plotspecs["saveas"] + "." + format, bbox_inches="tight")
                    os.chdir(cwd)



def plot_replica_maps_grid(dataset, plotspecs):
    """Plot contact maps for all replicas of each b-value in one grid"""
    cwd = os.getcwd()
    grid_dims = plotspecs["grid_dims"]
    for t in range(len(dataset.topologies)):
        names = dataset.top_names[t]
        for n in range(len(names)):
            # Plot whatever for a protein
            pairs = dataset.pairs[t][n]
            N = dataset.prot_sizes[t][n]
            print dataset.top_names[t][n]
            for j in range(len(dataset.b_values)):
                print "  b-values:", dataset.b_values[j]
                fig, axes = plt.subplots(*grid_dims, sharex=True, sharey=True, figsize=(12,10))
                if len(dataset.ydata[t][n][j]) > 0:
                    for rep in range(len(dataset.ydata[t][n][j])):
                        ax = axes[rep / grid_dims[0], rep % grid_dims[0]]

                        vals = dataset.ydata[t][n][j][0]
                        C = np.zeros((N, N))
                        for m in range(len(pairs)):
                            if m < dataset.prot_n_native[t][n]:
                                C[pairs[m, 1], pairs[m, 0]] = vals[m]
                            else:
                                C[pairs[m, 1], pairs[m, 0]] = -vals[m]

                        # plot native and non-native contacts in different colors
                        vmin, vmax = plotspecs["vminmax"]
                        pa = ax.pcolormesh(np.ma.array(C, mask=(C == 0)), cmap="bwr_r", vmin=vmin, vmax=vmax)

                        ax.annotate("rep = " + str(rep + 1),
                            xy=(0,0), xytext=plotspecs["xytext"],
                            bbox={"boxstyle":"square","facecolor":"w","edgecolor":"k"},
                            xycoords="axes fraction", textcoords="axes fraction")
                        ax.plot(np.arange(0, N), np.arange(0, N), 'k', lw=2)

                        ax.set_xlim(0, N)
                        ax.set_ylim(0, N)
                        ax.set_aspect("equal")

                    plt.subplots_adjust(wspace=0, hspace=0)
                    big_ax = fig.add_subplot(111)
                    big_ax.grid(False)
                    big_ax.set_axis_bgcolor('none')
                    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
                    big_ax.set_ylabel(plotspecs["ylabel"])
                    big_ax.set_xlabel(plotspecs["xlabel"])
                    big_ax.set_title(plotspecs["title"] + " b = " + dataset.b_values[j])

                    if not (plotspecs["saveas"] is None):
                        savedir = "{}/{}/b_{}/plots".format(dataset.topologies[t], 
                                dataset.top_names[t][n], dataset.b_values[j])
                        
                        if not os.path.exists(savedir):
                            os.mkdir(savedir)
                        os.chdir(savedir)
                        for format in plotspecs["saveas_formats"]:
                            plt.savefig(plotspecs["saveas"] + "." + format, bbox_inches="tight")
                        os.chdir(cwd)

def plot_bvalue_maps_grid(dataset, plotspecs):
    """Plot contact maps for all b-value in one grid. E.x. average contact map"""
    grid_dims = plotspecs["grid_dims"]
    for t in range(len(dataset.topologies)):
        names = dataset.top_names[t]
        for n in range(len(names)):
            # Plot whatever for a protein
            pairs = dataset.pairs[t][n]
            N = dataset.prot_sizes[t][n]
            print dataset.top_names[t][n]
            fig, axes = plt.subplots(*grid_dims, sharex=True, sharey=True, figsize=(12,10))
            for j in range(len(dataset.b_values)):
                print "  b-values:", dataset.b_values[j]
                ax = axes[j / grid_dims[0], j % grid_dims[0]]
                if len(dataset.ydata[t][n][j]) > 0:
                    vals = dataset.ydata[t][n][j][0]

                    C = np.zeros((N, N))
                    for m in range(len(pairs)):
                        if m < dataset.prot_n_native[t][n]:
                            C[pairs[m, 1], pairs[m, 0]] = vals[m]
                        else:
                            C[pairs[m, 1], pairs[m, 0]] = -vals[m]

                    # plot native and non-native contacts in different colors
                    vmin, vmax = plotspecs["vminmax"]
                    pa = ax.pcolormesh(np.ma.array(C, mask=(C == 0)), cmap="bwr_r", vmin=vmin, vmax=vmax)

                    ax.plot(np.arange(0, N + 1), np.arange(0, N + 1), 'k', lw=2)

                    ax.set_xlim(0, N)
                    ax.set_ylim(0, N)

                ax.annotate("b = " + dataset.b_values[j],
                    xy=(0,0), xytext=plotspecs["xytext"],
                    bbox={"boxstyle":"square","facecolor":"w","edgecolor":"k"},
                    xycoords="axes fraction", textcoords="axes fraction")
                ax.plot(np.arange(0, N), np.arange(0, N), 'k', lw=2)

            plt.subplots_adjust(wspace=0, hspace=0)
            big_ax = fig.add_subplot(111)
            big_ax.grid(False)
            big_ax.set_axis_bgcolor('none')
            big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
            big_ax.set_ylabel(plotspecs["ylabel"])
            big_ax.set_xlabel(plotspecs["xlabel"])
            big_ax.set_title(plotspecs["title"] + " " + dataset.top_names[t][n])
            #plt.colorbar(pa, big_ax, )

            if not (plotspecs["saveas"] is None):
                if not os.path.exists("plots"):
                    os.mkdir("plots")
                os.chdir("plots")
                for format in plotspecs["saveas_formats"]:
                    plt.savefig(plotspecs["saveas"].format(dataset.top_names[t][n]) + "." + format, bbox_inches="tight")
                os.chdir("..")






