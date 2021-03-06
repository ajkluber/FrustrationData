import os
import numpy as np
import matplotlib.pyplot as plt

import project_plotter
import project_util

import plotter.pmfutil

class dFdagg_Dataset(project_util.Dataset):
    def __init__(self, *args, **kwargs):
        project_util.Dataset.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        T_used = project_util.get_T_used()
        if os.path.exists(self.datapath + "T_{:.2f}_mid_bin.dat".format(T_used)):
            mid_bin = np.loadtxt(self.datapath + "T_{:.2f}_mid_bin.dat".format(T_used))
            FvsQ = np.loadtxt(self.datapath + "T_{:.2f}_F.dat".format(T_used))
        elif os.path.exists(self.datapath + "T_{:.1f}_mid_bin.dat".format(T_used)):
            mid_bin = np.loadtxt(self.datapath + "T_{:.1f}_mid_bin.dat".format(T_used))
            FvsQ = np.loadtxt(self.datapath + "T_{:.1f}_F.dat".format(T_used))

        xmin, xmax, F = plotter.pmfutil.find_U_TS_N(mid_bin, FvsQ)
        if len(xmin) == 0: 
            data = np.nan
        else:
            data = F(xmax[0]) - F(xmin[0])
        return data

if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]
    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25",
                "1.35", "1.45"]
    float_b = [ float(x) for x in b_values ]

    replicas = range(1,11)

    plotstyle = project_plotter.global_plot_style()

    mfpt_ylog_norm_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$\ln{\left(\frac{t_{mfpt}}{t_{mfpt}^{b=0}}\right)}$", "legend_loc":2,
            "ylim":(-1, 2)}
            #"ylim":(-1.5, 1.5)}
    mfpt_ylog_norm_plotspecs.update(plotstyle)

    mfpt_norm_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$\frac{t_{mfpt}}{t_{mfpt}^{b=0}}$", "legend_loc":2}
    mfpt_norm_plotspecs.update(plotstyle)

    msm_mfpt_ylog_norm_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$\ln{\left(\frac{t_{mfpt}}{t_{mfpt}^{b=0}}\right)}$", 
            "ylim":(-1, 2)}
            #"ylim":(-1.5, 1.5)}
    msm_mfpt_ylog_norm_plotspecs.update(plotstyle)
    msm_mfpt_ylog_norm_plotspecs["linestyle"] = "--"

    msm_mfpt_norm_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$\frac{t_{mfpt}}{t_{mfpt}^{b=0}}$"}
    msm_mfpt_norm_plotspecs.update(plotstyle)
    msm_mfpt_norm_plotspecs["linestyle"] = "--"

    prefactor_norm_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$\frac{\tau_{0}}{\tau_{0}^{b=0}}$", "legend_loc":2}
    prefactor_norm_plotspecs.update(plotstyle)

    tp_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$t_{tp}$ (frames)"}
    tp_plotspecs.update(plotstyle)

    tp_norm_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$\frac{t_{tp}}{t_{tp}^{b=0}}$"}
    tp_norm_plotspecs.update(plotstyle)

    #saveas = "two_colm_mfpt_tpt"
    saveas = "two_colm_pref_tpt"

    # gather data
    msm_mfpt_data = project_util.Dataset(topologies, top_names, b_values, replicas, "msm/tfold_pcca.dat")
    mfpt_data = project_util.Dataset(topologies, top_names, b_values, replicas, "Qtanh_0_05_folding_time/folding_mean")
    tp_data = project_util.Dataset(topologies, top_names, b_values, replicas, "Qtanh_0_05_transit_time/forward_mean")
    dFdagg_raw_data = dFdagg_Dataset(topologies, top_names, b_values, replicas, "Qtanh_0_05_profile/", dont_avg=True)

    mfpt_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            mfpt_0 = mfpt_data.data[t][n][0]
            avg_mfpt_0 = np.mean(mfpt_0[~np.isnan(mfpt_0)])
            for b in range(len(b_values)):
                mfpt_b = mfpt_data.data[t][n][b]
                for rep in range(len(mfpt_b)):
                    if not np.isnan(mfpt_b[rep]):
                        mfpt_norm_data.data[t][n][b][rep] = mfpt_b[rep]/avg_mfpt_0
    mfpt_norm_data._calc_repavg()

    # normalized logarithm
    mfpt_ylog_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            mfpt_0 = mfpt_data.data[t][n][0]
            avg_mfpt_0 = np.mean(mfpt_0[~np.isnan(mfpt_0)])
            for b in range(len(b_values)):
                mfpt_b = mfpt_data.data[t][n][b]
                for rep in range(len(mfpt_b)):
                    if not np.isnan(mfpt_b[rep]):
                        mfpt_ylog_norm_data.data[t][n][b][rep] = np.log(mfpt_b[rep]/avg_mfpt_0)
    mfpt_ylog_norm_data._calc_repavg()

    msm_mfpt_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            mfpt_0 = msm_mfpt_data.data[t][n][0]
            avg_mfpt_0 = np.mean(mfpt_0[~np.isnan(mfpt_0)])
            for b in range(len(b_values)):
                mfpt_b = msm_mfpt_data.data[t][n][b]
                for rep in range(len(mfpt_b)):
                    if not np.isnan(mfpt_b[rep]):
                        msm_mfpt_norm_data.data[t][n][b][rep] = mfpt_b[rep]/avg_mfpt_0
    msm_mfpt_norm_data._calc_repavg()

    msm_mfpt_ylog_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            mfpt_0 = msm_mfpt_data.data[t][n][0]
            avg_mfpt_0 = np.mean(mfpt_0[~np.isnan(mfpt_0)])
            for b in range(len(b_values)):
                mfpt_b = msm_mfpt_data.data[t][n][b]
                for rep in range(len(mfpt_b)):
                    if not np.isnan(mfpt_b[rep]):
                        msm_mfpt_ylog_norm_data.data[t][n][b][rep] = np.log(mfpt_b[rep]/avg_mfpt_0)
    msm_mfpt_ylog_norm_data._calc_repavg()

    prefactor_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                mfpt_b = mfpt_data.data[t][n][b]
                dF_b = dFdagg_raw_data.data[t][n][b]
                for rep in range(len(mfpt_b)):
                    if (not np.isnan(mfpt_b[rep])) and (not np.isnan(dF_b[rep])):
                        prefactor_data.data[t][n][b][rep] = mfpt_b[rep]/np.exp(dF_b[rep])
    prefactor_data._calc_repavg()

    prefactor_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                pre_0 = prefactor_data.data[t][n][0]
                pre_unfrust = np.mean(pre_0[~np.isnan(pre_0)])
                pre_b = prefactor_data.data[t][n][b]
                for rep in range(len(pre_b)):
                    if not np.isnan(pre_b[rep]):
                        prefactor_norm_data.data[t][n][b][rep] = pre_b[rep]/pre_unfrust
    prefactor_norm_data._calc_repavg()

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
    tp_norm_data._calc_repavg()

    semilog = False

    #all_datasets = [[mfpt_norm_data, msm_mfpt_norm_data], [tp_norm_data]]
    #all_plotspecs = [[mfpt_norm_plotspecs, msm_mfpt_norm_plotspecs], [tp_norm_plotspecs]]
    all_datasets = [[prefactor_norm_data], [tp_norm_data]]
    all_plotspecs = [[prefactor_norm_plotspecs], [tp_norm_plotspecs]]

    panel_label = ["A", "B"]

    # plot two colm panel
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))
    for d in [0,1]: 
        ax = axes[d]
        for z in range(len(all_datasets[d])): 
            dataset = all_datasets[d][z]
            plotspecs = all_plotspecs[d][z]
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
                            rep_data = dataset.data[t][n][j]
                            good = ~np.isnan(rep_data)
                            avg_val = np.mean(rep_data[good])
                            std_val = np.std(rep_data[good])/np.sqrt(float(np.sum(good)))

                        # plot average value with error bars 
                        if not np.isnan(avg_val):
                            b_line.append(b)
                            avg_val_line.append(avg_val)
                            std_dev_line.append(std_val)

                    if plotspecs.has_key("legend_key"):
                        ax.errorbar(b_line, avg_val_line, yerr=std_dev_line,
                                color=plotspecs["color"][t][n],
                                marker=plotspecs["markers"][t], ls=plotspecs["linestyle"],
                                label=plotspecs["legend_key"][t][n], lw=2, ms=8)
                    else:
                        ax.errorbar(b_line, avg_val_line, yerr=std_dev_line,
                                color=plotspecs["color"][t][n], ls=plotspecs["linestyle"],
                                marker=plotspecs["markers"][t], lw=2, ms=8)

            if semilog:
                if d in [0,1]:
                    ax.semilogy()
            if z == 0:
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')

                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(14) 
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(14) 

                if plotspecs.has_key("ylim"):
                    ax.set_ylim(plotspecs["ylim"])

                if plotspecs.has_key("legend_loc"):
                    ax.legend(loc=plotspecs["legend_loc"], fontsize=16)

                #ax.set_ylabel(plotspecs["ylabel"], fontsize=20, rotation="horizontal")
                ax.annotate(plotspecs["ylabel"], fontsize=28, xy=(0,0), xytext=(-.224, 0.5),
                        xycoords="axes fraction", textcoords="axes fraction")

                ax.set_xlabel(plotspecs["xlabel"], fontsize=20)
                plt.subplots_adjust(right=0.98, left=0.02, wspace=0.25)
            
                ax.annotate(panel_label[d], xy=(0,0), xytext=(-0.19, 1.05), fontsize=30,
                    xycoords="axes fraction", textcoords="axes fraction")

    #plt.show()

    for format in plotspecs["saveas_formats"]:
        fig.savefig("plots/{}".format(saveas) + "." + format, bbox_inches="tight")
    plt.show()

