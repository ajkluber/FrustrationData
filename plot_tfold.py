import os
import numpy as np
import matplotlib.pyplot as plt

import project_plotter
import project_util

import plotter.pmfutil

class Dataset(project_util.Dataset):
    def __init__(self, *args, **kwargs):
        project_util.Dataset.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        with open(self.datapath, "r") as fin:
            data = float(fin.read())
        return data


class dFstab_Dataset(project_util.Dataset):
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
        if len(xmin) < 2: 
            data = np.nan
        else:
            data = F(xmin[1]) - F(xmin[0])
        return data

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

    replicas = range(1, 11)

    plotstyle = project_plotter.global_plot_style()

    coordname = "Qtanh_0_05"
    datapath = coordname + "_folding_time/folding_mean"

    mfpt_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$t_{mfpt}$ (frames)", 
            "title":"MFPT along Q", "legend_loc":2,
            "saveas":"repavg_mfpt_1D"}
    mfpt_plotspecs.update(plotstyle)

    mfpt_ylog_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":r"$\log\frac{t_{mfpt}}{t_{unit}}$",
            "title":"MFPT along Q", "legend_loc":6, "ylims":(5.5,10), 
            "saveas":"repavg_mfpt_1D_ylog"}
    mfpt_ylog_plotspecs.update(plotstyle)

    mfpt_ylog_norm_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":r"$\log\frac{t_{mfpt}}{t^{b=0}_{mfpt}}$", 
            "title":"MFPT along Q", "legend_loc":2,
            "saveas":"repavg_mfpt_1D_ylog_norm"}
    mfpt_ylog_norm_plotspecs.update(plotstyle)

    mfpt_norm_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":r"$t_{mfpt}$ / $t_{mfpt}^{b=0}$", 
            "title":"MFPT along Q", "legend_loc":6, "ylims":(0, 10),
            "saveas":"repavg_mfpt_1D_norm"}
    mfpt_norm_plotspecs.update(plotstyle)

    dFstab_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":r"$\Delta F^{\circ}$ $(k_B T)$", 
            "title":"Stability", "legend_loc":3,
            "saveas":"repavg_dFstab_1D"}
    dFstab_plotspecs.update(plotstyle)

    dFdagg_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":r"$\Delta F^{\dagger}$ $(k_B T)$", 
            "title":"Barrier height", "legend_loc":3, "ylims":(0,5),
            "saveas":"repavg_dFdagg_1D"}
    dFdagg_plotspecs.update(plotstyle)

    prefactor_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":r"$t_0$ (frames)", 
            "title":"Folding prefactor","legend_loc":2, "ylims":(0, 900),
            "saveas":"repavg_prefactor_1D"}
    prefactor_plotspecs.update(plotstyle)

    prefactor_norm_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":r"$t_0$ / $t_0^{b=0}$", 
            "title":"Folding prefactor", "legend_loc":2, "ylims":(0, 5),
            "saveas":"repavg_prefactor_1D_norm"}
    prefactor_norm_plotspecs.update(plotstyle)

    # gather data
    mfpt_data = project_util.Dataset(topologies, top_names, b_values, replicas, datapath)
    dFstab_raw_data = dFstab_Dataset(topologies, top_names, b_values, replicas, coordname + "_profile/", dont_avg=True)
    dFdagg_raw_data = dFdagg_Dataset(topologies, top_names, b_values, replicas, coordname + "_profile/", dont_avg=True)

    #project_util.save_avg_plot_data(mfpt_data, "mfpt_Q")
    #project_util.save_avg_plot_data(dFdagg_raw_data, "dFdagg")
    #raise SystemExit

    # normalize data 
    mfpt_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            mfpt_0 = np.array(mfpt_data.data[t][n][0])
            avg_mfpt_0 = np.mean(mfpt_0[~np.isnan(mfpt_0)])
            for b in range(len(b_values)):
                mfpt_b = mfpt_data.data[t][n][b]
                for rep in range(len(mfpt_b)):
                    if not np.isnan(mfpt_b[rep]):
                        mfpt_norm_data.data[t][n][b][rep] = mfpt_b[rep]/avg_mfpt_0
    mfpt_norm_data._calc_repavg()

    mfpt_ylog_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            mfpt_0 = np.array(mfpt_data.data[t][n][0])
            avg_mfpt_0 = np.mean(mfpt_0[~np.isnan(mfpt_0)])
            for b in range(len(b_values)):
                mfpt_b = mfpt_data.data[t][n][b]
                for rep in range(len(mfpt_b)):
                    if not np.isnan(mfpt_b[rep]):
                        mfpt_ylog_norm_data.data[t][n][b][rep] = np.log(mfpt_b[rep]/avg_mfpt_0)
    mfpt_ylog_norm_data._calc_repavg()

    mfpt_ylog_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                mfpt_ylog_b = mfpt_data.data[t][n][b]
                for rep in range(len(mfpt_ylog_b)):
                    if not np.isnan(mfpt_ylog_b[rep]):
                        mfpt_ylog_data.data[t][n][b][rep] = np.log(mfpt_ylog_b[rep])


    # deconvolute effects of barrier height
    prefactor_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                mfpt_b = mfpt_data.data[t][n][b]
                dF_b = dFdagg_raw_data.data[t][n][b]
                for rep in range(len(mfpt_b)):
                    if (not np.isnan(mfpt_b[rep])) and (not np.isnan(dF_b[rep])):
                        prefactor_data.data[t][n][b][rep] = mfpt_b[rep]/np.exp(dF_b[rep])

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

    dFstab_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                dF_b = dFstab_raw_data.data[t][n][b]
                for rep in range(len(dF_b)):
                    if not np.isnan(dF_b[rep]):
                        dFstab_data.data[t][n][b][rep] = dF_b[rep]
    dFstab_data._calc_repavg()

    dFdagg_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                dF_b = dFdagg_raw_data.data[t][n][b]
                for rep in range(len(dF_b)):
                    if not np.isnan(dF_b[rep]):
                        dFdagg_data.data[t][n][b][rep] = dF_b[rep]
    dFdagg_data._calc_repavg()


    # plot data
    project_plotter.plot_data(mfpt_ylog_norm_data, mfpt_ylog_plotspecs)
    project_plotter.plot_data(mfpt_ylog_data, mfpt_ylog_plotspecs)
    project_plotter.plot_data(mfpt_data, mfpt_plotspecs)
    project_plotter.plot_data(mfpt_norm_data, mfpt_norm_plotspecs)
    project_plotter.plot_data(dFstab_data, dFstab_plotspecs)
    project_plotter.plot_data(dFdagg_data, dFdagg_plotspecs)
    project_plotter.plot_data(prefactor_data, prefactor_plotspecs)
    project_plotter.plot_data(prefactor_norm_data, prefactor_norm_plotspecs)
    plt.show()

