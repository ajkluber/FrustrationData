import os
import numpy as np
import matplotlib.pyplot as plt

import project_util, project_plotter

if __name__ == "__main__":
    #topologies = ["alpha"]
    #top_names = [["1imq"]]

    #topologies = ["beta", "alpha"]
    #top_names = [["2akk"], ["1imq"]]

    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    replicas = range(1,11)

    plotstyle = project_plotter.global_plot_style()

    Tf_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"Temperature", 
            "legend_loc":3}
    Tf_plotspecs.update(plotstyle)

    Tg_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"Temperature", 
            "legend_loc":2}
    Tg_plotspecs.update(plotstyle)

    Tg_over_Tf_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":r"$\frac{T_g}{T_f}$", 
            "legend_loc":2, "saveas":None}
    Tg_over_Tf_plotspecs.update(plotstyle)

    Tf_data = project_util.Dataset(topologies, top_names, b_values, replicas, "Qtanh_0_05_profile/T_used.dat")
    Tg_data = project_util.Dataset(topologies, top_names, b_values, replicas, "Tg_calc/Tg_Enonnative.dat")

    Tg_over_Tf_data = project_util.Dataset(topologies, top_names, b_values, replicas, "", empty=True)
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            for b in range(len(b_values)):
                Tf_b = Tf_data.data[t][n][b]
                Tg_b = Tg_data.data[t][n][b]
                for rep in range(len(Tf_b)):
                    if not np.isnan(Tf_b[rep]) and not np.isnan(Tg_b[rep]):
                        Tg_over_Tf_data.data[t][n][b][rep] = Tg_b[rep]/Tf_b[rep]
    Tg_over_Tf_data._calc_repavg()

    project_plotter.plot_data(Tg_over_Tf_data, Tg_over_Tf_plotspecs)
    plt.show()

    raise SystemExit
    Tf_data = get_Tf_data(topologies, top_names, b_values, replicas, "Qtanh_0_05")
    Tg_data = get_Tg_data(topologies, top_names, b_values, replicas)
    plt.figure()
    # plot disorder-average temperatures
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            first_time = True
            avg_Tf_unfrust = np.mean(Tf_all[t][n][0])
            for b in range(len(float_b)):
                if len(Tf_all[t][n][b]) > 0:
                    # plot disorder average
                    avg_Tf = np.mean(Tf_all[t][n][b])
                    avg_Tg = np.mean(Tg_all[t][n][b])
                    std_Tf = np.std(Tf_all[t][n][b])
                    std_Tg = np.std(Tg_all[t][n][b])
                    if first_time:
                        #plt.plot(float_b[b], avg_Tf/avg_Tf_unfrust, marker=symbols[t],
                        #        color=color_codes[t][n], 
                        #        label=top_names[t][n] + " (" + str(prot_sizes[t][n]) + ")")
                        plt.errorbar(float_b[b], avg_Tf, yerr=std_Tf, marker=symbols[t],
                                color=color_codes[t][n], 
                                label=top_names[t][n] + " (" + str(prot_sizes[t][n]) + ")")
                        plt.errorbar(float_b[b], avg_Tg, yerr=std_Tg, marker=symbols[t],
                                color=color_codes[t][n])

                        #plt.plot(float_b[b], avg_Tg, marker=symbols[t],
                        #        color=color_codes[t][n])
                        first_time = False
                    else:
                        #plt.plot(float_b[b], avg_Tf/avg_Tf_unfrust, marker=symbols[t], color=color_codes[t][n])
                        #plt.plot(float_b[b], avg_Tg, marker=symbols[t], color=color_codes[t][n])
                        plt.errorbar(float_b[b], avg_Tf, yerr=std_Tf, marker=symbols[t], color=color_codes[t][n])
                        plt.errorbar(float_b[b], avg_Tg, yerr=std_Tg, marker=symbols[t], color=color_codes[t][n])

    plt.legend(loc=(0.05, 0.25), fontsize=16)
    plt.xlabel("Frustration $b$")
    plt.ylabel("Temperature")
    ymin, ymax = plt.ylim()
    plt.ylim(0, ymax)
    plt.savefig("plots/Tf_and_Tg.png")
    plt.savefig("plots/Tf_and_Tg.pdf")

    plt.figure()
    # plot disorder-average temperatures
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            first_time = True
            avg_Tf_unfrust = np.mean(Tf_all[t][n][0])
            for b in range(len(float_b)):
                if len(Tf_all[t][n][b]) > 0:
                    # plot disorder average
                    avg_Tf = np.mean(Tf_all[t][n][b])
                    avg_Tg = np.mean(Tg_all[t][n][b])
                    if first_time:
                        plt.plot(float_b[b], avg_Tg/avg_Tf, marker=symbols[t], color=color_codes[t][n],
                                label=top_names[t][n] + " (" + str(prot_sizes[t][n]) + ")")
                        first_time = False
                    else:
                        plt.plot(float_b[b], avg_Tg/avg_Tf, marker=symbols[t], color=color_codes[t][n])

    plt.legend(loc=2, fontsize=16)
    plt.xlabel("Frustration $b$")
    plt.ylabel("Glassiness $\\left(\\frac{T_g}{T_f}\\right)$")
    plt.savefig("plots/Tf_over_Tg.pdf")
    plt.ylim(0, 1)
    plt.show()
