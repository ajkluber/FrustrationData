import os
import glob
import numpy as np
import matplotlib.pyplot as plt

import project_util
import project_plotter

class Dataset(project_util.Dataset):
    def __init__(self, *args, **kwargs):
        self.unfolded_Q = None
        self.barrier_Q = None
        self.datapath_x = kwargs.pop("datapath_x")
        project_util.Dataset.__init__(self, *args, **kwargs)

    def get_replica_data(self):
        if self.unfolded_Q is None:
            # determine which Q bin to use
            U = np.loadtxt("Qtanh_0_05_profile/minima.dat")[0]
            self.unfolded_Q = U
            #TS = np.loadtxt("Qtanh_0_05_profile/maxima.dat")[0]
            #self.barrier_Q = TS
        else:
            U = self.unfolded_Q
            #TS = self.barrier_Q

        Q_mid_bin = np.load(self.datapath_x)
        Q_bin_idx = np.argmin((Q_mid_bin - U)**2)

        Rg_vs_Q = np.load(self.datapath)
        data = Rg_vs_Q[Q_bin_idx]

        return data

def plot_Rg_scaled_by_powerlaw():
    
    # Plot radius of gyration normalized by length
    Rg_scaled_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$R_g$ / $N^{1.1}$ (nm)", 
            "ylims":(0.005,0.015), "title":"Scaled Radius of Gyration", "legend_loc":3,
            "saveas":"Rg_U_scaled_N_vs_b"}
    Rg_scaled_plotspecs.update(plotstyle)

    Rg_scaled_data = project_util.Dataset(topologies, top_names, b_values, replicas,"", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                Rg_b = Rg_data.data[t][n][b]
                for rep in range(len(Rg_b)):
                    if not np.isnan(Rg_b[rep]):
                        #Rg_scaled_data.data[t][n][b].append(Rg_b[rep]/(float(prot_sizes[t][n])))
                        Rg_scaled_data.data[t][n][b].append(Rg_b[rep]/(float(prot_sizes[t][n])**1.1)) # My weird scaling


    # Plot radius of gyration scaled by length for a self-avoiding walk SAW. N^(3/5)
    Rg_SAW_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":"$R_g$ / $N^{3/5}$ (nm)", "ylims":(0.0, 0.14),
            "title":"SAW scaled Radius of gyration", "legend_loc":3,
            "saveas":"Rg_U_SAW_vs_b"}
    Rg_SAW_plotspecs.update(plotstyle)

    Rg_SAW_data = project_util.Dataset(topologies, top_names, b_values, replicas,"", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                Rg_b = Rg_data.data[t][n][b]
                for rep in range(len(Rg_b)):
                    if not np.isnan(Rg_b[rep]):
                        Rg_SAW_data.data[t][n][b].append(Rg_b[rep]/(float(prot_sizes[t][n])**(3./5))) # SAW

    project_plotter.plot_data(Rg_scaled_data, Rg_scaled_plotspecs)
    project_plotter.plot_data(Rg_SAW_data, Rg_SAW_plotspecs)


if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    N = np.array([63., 86., 58., 74., 48.])

    plotstyle = project_plotter.global_plot_style()

    coordname = "Qtanh_0_05"
    datapath_x = "Rg_vs_" + coordname + "/mid_bin.npy"
    datapath = "Rg_vs_" + coordname + "/Rg_vs_bin.npy"


    replicas = range(1,11)
    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25", "1.35", "1.45"]

    Rg_data = Dataset(topologies, top_names, b_values, replicas, datapath,
            datapath_x=datapath_x)


    # Plot radius of gyration for unfolded state
    Rg_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$R_g$ (nm)", "ylims": (0,2),
            "title":"Radius of Gyration", "legend_loc":3, "saveas":"Rg_U_vs_b"}
    Rg_plotspecs.update(plotstyle)


    # Plot radius of gyration normalized to unfrustrated case
    Rg_norm_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"$R_g$ / $R_g^{b=0}$", 
            "ylims":(0.5, 1.04), "legend_loc":3, "saveas":"Rg_U_norm_vs_b"}
    Rg_norm_plotspecs.update(plotstyle)


    Rgmin_over_Rg_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$\frac{R_g^{min}}{R_g}$", "ylims":(0.3, 1), 
            "legend_loc":4, "saveas":"Rgmin_over_Rg_vs_b"}
    Rgmin_over_Rg_plotspecs.update(plotstyle)

    eta_pack_plotspecs = {"xlabel":"Frustration ($b$)", 
            "ylabel":r"$\eta$", "ylims":(0, 1), 
            "legend_loc":2, "saveas":"eta_pack_vs_b"}
    eta_pack_plotspecs.update(plotstyle)

    Rg_norm_data = project_util.Dataset(topologies, top_names, b_values, replicas,"", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            Rg_0 = Rg_data.data[t][n][0]
            avg_Rg_0 = np.mean(Rg_0[~np.isnan(Rg_0)])
            for b in range(len(b_values)):
                Rg_b = Rg_data.data[t][n][b]
                for rep in range(len(Rg_b)):
                    if not np.isnan(Rg_b[rep]):
                        Rg_norm_data.data[t][n][b][rep] = Rg_b[rep]/avg_Rg_0

    Rg_min = lambda N: np.sqrt(3./5)*(0.35*((1. - 0.3)*N - 1)**(1./3)) # Rg of solid sphere of N monomers.

    coeff1, coeff2 = -4.26638807901, 1.10395088206
    Rg_coil = lambda N: np.exp(coeff1)*(N**coeff2) # Rg of coil 


    Rgmin_over_Rg = project_util.Dataset(topologies, top_names, b_values, replicas,"", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            Rgmin = Rg_min(float(prot_sizes[t][n]))
            for b in range(len(b_values)):
                Rg_b = Rg_data.data[t][n][b]
                for rep in range(len(Rg_b)):
                    if not np.isnan(Rg_b[rep]):
                        Rgmin_over_Rg.data[t][n][b][rep] = Rgmin/Rg_b[rep]
    Rgmin_over_Rg._calc_repavg()

    eta_pack = project_util.Dataset(topologies, top_names, b_values, replicas,"", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            Rgmin = Rg_min(float(prot_sizes[t][n]))
            #Rgcoil = Rg_coil(float(prot_sizes[t][n]))
            Rgcoil = np.max(Rg_data.data[t][n][0])
            for b in range(len(b_values)):
                Rg_b = Rg_data.data[t][n][b]
                for rep in range(len(Rg_b)):
                    if not np.isnan(Rg_b[rep]):
                        eta_pack.data[t][n][b][rep] = (Rgcoil - Rg_b[rep])/(Rgcoil - Rgmin)
    eta_pack._calc_repavg()

    # correct error bars on packing fraction
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            Rgmin = Rg_min(float(prot_sizes[t][n]))
            temp_std = []
            for b in range(len(b_values)):
                Rg_b = Rg_data.data[t][n][b]
                good = ~np.isnan(Rg_b)
                if np.sum(good) > 1:
                    Rg_avg = np.mean(Rg_b[good])
                    Rg_SE = np.std(Rg_b[good])/np.sqrt(len(Rg_b[good]))
                    Rg_percent_error = Rg_SE/Rg_avg
                    new_error = Rg_percent_error*eta_pack.avgdata[t][n][b]
                    temp_std.append(new_error)
                else:
                    temp_std.append(np.nan)
            eta_pack.stddata[t][n] = np.array(temp_std)

    # correct error bars on the reciprocal
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            Rgmin = Rg_min(float(prot_sizes[t][n]))
            temp_std = []
            for b in range(len(b_values)):
                Rg_b = Rg_data.data[t][n][b]
                good = ~np.isnan(Rg_b)
                if np.sum(good) > 1:
                    Rg_avg = np.mean(Rg_b[good])
                    Rg_SE = np.std(Rg_b[good])/np.sqrt(len(Rg_b[good]))
                    Rg_percent_error = Rg_SE/Rg_avg
                    new_error = (Rgmin/Rg_avg)*Rg_percent_error
                    temp_std.append(new_error)
                else:
                    temp_std.append(np.nan)
            Rgmin_over_Rg.stddata[t][n] = np.array(temp_std)

    #project_plotter.plot_data(Rg_data, Rg_plotspecs)
    #project_plotter.plot_data(Rg_norm_data, Rg_norm_plotspecs)
    #project_plotter.plot_data(Rgmin_over_Rg, Rgmin_over_Rg_plotspecs)
    project_plotter.plot_data(eta_pack, eta_pack_plotspecs)
    plt.show()

    project_util.save_avg_plot_data(Rg_data, "Rg_raw")
    project_util.save_avg_plot_data(Rgmin_over_Rg, "Rg_min_over_Rg")

