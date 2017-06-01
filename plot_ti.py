import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

import pyemma.msm as msm

def get_timescale_data(topologies, top_names, b_values, replicas, coordname, n_timescales):
    msm_ti = []
    for t in range(len(topologies)):
        os.chdir(topologies[t])
        names = top_names[t]
        msm_ti_top = []
        for name in names:
            msm_ti_name = []
            print name
            os.chdir(name)
            for b in b_values:
                print "b-value: ", b

                msm_ti_rep = []
                for rep in replicas:
                    os.chdir("b_{}/replica_{}".format(b, rep))
                    #if os.path.exists("msm/timescales.npy"):
                    #    # load in implied timescales
                    #    lagtimes = np.load("msm/lagtimes.npy")
                    #    its_vs_lag = np.load("msm/timescales.npy")
                    #    msm_ti_rep.append(its_vs_lag)
                    #elif os.path.exists("msm/msm.pkl"):
                    if True:
                        with open("msm/msm.pkl", "rb") as fhandle:
                            msm_info = pickle.load(fhandle)

                        lagtimes = msm_info["lagtimes"]
                        its_vs_lag = np.zeros((len(lagtimes), n_timescales))
                        for i in range(len(lagtimes)):
                            lag = lagtimes[i]
                            model = msm.markov_model(msm_info[str(lag)])
                            its_vs_lag[i,:] = lag*model.timescales()[:n_timescales]

                        # save implied timescale info
                        np.save("msm/lagtimes.npy", lagtimes)
                        np.save("msm/timescales.npy", its_vs_lag)
                        msm_ti_rep.append(its_vs_lag)

                    os.chdir("../..")
                    if b == "0.00":
                        # if unfrustrated go-model only run one replica.
                        break
                msm_ti_name.append(msm_ti_rep)
            os.chdir("..")
            msm_ti_top.append(msm_ti_name)
        msm_ti.append(msm_ti_top)
        os.chdir("..")
    return msm_ti

def plot_protein_grid_ti_noavg(msm_ti, topologies, top_names, float_b, replicas, lag_idx, time_idx, semilogy, ti_colors):
    fig, axes = plt.subplots(1, 5, figsize=(20, 4), sharex=True, sharey=True)
    counter = 0
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            lag_idx = lag_idxs[t][n]
            ax = axes[counter]
            for b in range(len(float_b)):
                if len(msm_ti[t][n][b]) > 0:
                    for r in range(len(replicas)):
                        rep_ti = msm_ti[t][n][b][r]
                        if (b == 0) and (r == 0):
                            ax.annotate(top_names[t][n], xy=(0,0), xytext=(0.4,0.9),
                                    textcoords="axes fraction", xycoords="axes fraction", fontsize=16)

                        if semilogy:
                            ax.plot(float_b[b], np.log(rep_ti[lag_idx][time_idx]), ms=3, marker='o', color=ti_colors[b % len(ti_colors)])
                        else:
                            ax.plot(float_b[b], rep_ti[lag_idx][time_idx], ms=3, marker='o', color=ti_colors[b % len(ti_colors)])
            ax.set_xticklabels([ str(x) for x in ax.get_xticks() ], rotation=90)
            counter += 1

    fig.suptitle("MSM timescale $t_" + str(time_idx + 1) + "$", fontsize=18)

    big_ax = fig.add_subplot(111)
    big_ax.grid(False)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    big_ax.spines['top'].set_visible(False)
    big_ax.spines['bottom'].set_visible(False)
    big_ax.spines['right'].set_visible(False)
    big_ax.spines['left'].set_visible(False)

    if semilogy:
        ylabel = "$\\log t_" + str(time_idx + 1) + "$"
    else:
        ylabel = "$t_" + str(time_idx + 1) + "$"

    big_ax.annotate(ylabel, fontsize=18, rotation="vertical",
        xy=(0,0), xytext=(-0.08, 0.55),
        textcoords="axes fraction", xycoords="axes fraction")

    big_ax.annotate("Frustration (b)", fontsize=18,
        xy=(0,0), xytext=(0.45, -0.25),
        textcoords="axes fraction", xycoords="axes fraction")

    if semilogy:
        fig.savefig("grid_t" + str(time_idx + 1) + "_noavg_semilogy.png")
        fig.savefig("grid_t" + str(time_idx + 1) + "_noavg_semilogy.pdf")
    else:
        fig.savefig("grid_t" + str(time_idx + 1) + "_noavg.png")
        fig.savefig("grid_t" + str(time_idx + 1) + "_noavg.pdf")

def plot_ti_repavg(msm_ti, topologies, top_names, float_b, replicas, lag_idx, time_idx, semilogy, ti_colors, ylims=False):
    # plot disorder-average timescales.
    plt.figure()
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            first_time = True
            for b in range(len(float_b)):
                if len(msm_ti[t][n][b]) > 0:
                    b_ti = np.array([ msm_ti[t][n][b][r][lag_idx][time_idx] for r in range(len(msm_ti[t][n][b])) ])

                    if semilogy:
                        avg_ti = np.mean(np.log(b_ti))
                        err_ti = np.std(np.log(b_ti))/np.sqrt(float(len(b_ti)))
                    else:
                        avg_ti = np.mean(b_ti)
                        err_ti = np.std(b_ti)/np.sqrt(float(len(b_ti)))

                    if first_time:
                        plt.errorbar(float_b[b], avg_ti, yerr=err_ti, ms=4, marker=symbols[t], color=color_codes[t][n], label=top_names[t][n])
                        first_time = False
                    else:
                        plt.errorbar(float_b[b], avg_ti, yerr=err_ti, ms=4, marker=symbols[t], color=color_codes[t][n])

    if semilogy:
        ylabel = "$\\log \\frac{t_" + str(time_idx + 1) + "}{t_{unit}}$"
    else:
        ylabel = "$t_" + str(time_idx + 1) + "$"

    if ylims:
        plt.ylim(*ylims)

    plt.legend(loc=2)
    #if (time_idx == 0) and not semilogy:
    #    plt.legend(loc=6)
    #else:
    #    plt.legend(loc=4)

    plt.title("MSM timescale")
    plt.xlabel("Frustration $b$")
    plt.ylabel(ylabel)
    if semilogy:
        plt.savefig("repavg_msm_t" + str(time_idx + 1) + "_semilogy.png")
        plt.savefig("repavg_msm_t" + str(time_idx + 1) + "_semilogy.pdf")
    else:
        plt.savefig("repavg_msm_t" + str(time_idx + 1) + ".png")
        plt.savefig("repavg_msm_t" + str(time_idx + 1) + ".pdf")

def plot_ti_repavg_norm(msm_ti, topologies, top_names, float_b, replicas, lag_idx, time_idx, semilogy, ti_colors):

    # plot disorder-average timescales.
    plt.figure()
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            first_time = True
            
            for b in range(len(float_b)):
                unfrust_b_ti = np.array([ msm_ti[t][n][0][r][lag_idx][time_idx] for r in range(len(msm_ti[t][n][0])) ])
                if len(msm_ti[t][n][b]) > 0:
                    b_ti = np.array([ msm_ti[t][n][b][r][lag_idx][time_idx] for r in range(len(msm_ti[t][n][b])) ])

                    if semilogy:
                        avg_ti = np.mean(np.log(b_ti)) - np.mean(np.log(unfrust_b_ti))
                        #err_ti = np.std(np.log(b_ti))/np.sqrt(float(len(b_ti)))
                    else:
                        avg_ti = np.mean(b_ti)/np.mean(unfrust_b_ti)
                        #err_ti = np.std(b_ti)/np.sqrt(float(len(b_ti)))

                    if first_time:
                        #plt.errorbar(float_b[b], avg_ti, yerr=err_ti, ms=4, marker=symbols[t], color=color_codes[t][n], label=top_names[t][n])
                        plt.plot(float_b[b], avg_ti, ms=6, marker=symbols[t], color=color_codes[t][n], label=top_names[t][n])
                        first_time = False
                    else:
                        #plt.errorbar(float_b[b], avg_ti, yerr=err_ti, ms=4, marker=symbols[t], color=color_codes[t][n])
                        plt.errorbar(float_b[b], avg_ti, ms=6, marker=symbols[t], color=color_codes[t][n])

    if semilogy:
        ylabel = "$\\log \\frac{t_" + str(time_idx + 1) + "}{t^{b=0}_" + str(time_idx + 1) + "}$"
    else:
        ylabel = "$\\frac{t_" + str(time_idx + 1) + "}{t^{b=0}_" + str(time_idx + 1) + "}$"

    plt.legend(loc=2)
    plt.title("MSM timescale")
    plt.xlabel("Frustration $b$")
    plt.ylabel(ylabel)
    if semilogy:
        plt.savefig("repavg_msm_t" + str(time_idx + 1) + "_semilogy_norm.png")
        plt.savefig("repavg_msm_t" + str(time_idx + 1) + "_semilogy_norm.pdf")
    else:
        plt.savefig("repavg_msm_t" + str(time_idx + 1) + "_norm.png")
        plt.savefig("repavg_msm_t" + str(time_idx + 1) + "_norm.pdf")

if __name__ == "__main__":
    colorcycle = ["#5DA5DA", "#FAA43A", "#60BD68", "#F17CB0", "#B2912F", "#B276B2", "#DECF3F", "#F15854", "#4D4D4D"]

    #topologies = ["alpha", "beta"]
    #top_names = [["1imq"], ["2akk"]]

    #topologies = ["alpha", "mixed", "beta"]
    #symbols = ["^", "s", "o"]
    #top_names = [["1r69", "1imq"], ["1e0g"], ["1fmk", "2akk"]]
    #prot_sizes = [[63, 86], [48], [58, 74]]
    #color_codes = [["#5DA5DA", "#F17CB0"], ["#FAA43A"], ["#60BD68", "#B2912F"]]
    #lag_idxs = [[7, 9], [7], [8, 9]]

    #topologies = ["mixed"]
    #top_names = [["1e0g"]]

    #topologies = ["alpha"]
    #top_names = [["1r69"]]

    #topologies = ["beta"]
    #top_names = [["2akk"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]

    replicas = range(1,11)


    topologies = ["alpha"]
    top_names = [["1imq"]]
    replicas = [6]
    b_values = ["0.01"]

    coordname = "Qtanh_0_05"
    n_timescales = 10
    lagtime = 200
    lag_idx = 7
    time_idx = 1
    semilogy = False 
    ti_colors = ['r','g','b','salmon', 'darkgreen']
    lags = np.array([1,2,5,10,20,50,100,200,300,400,500,600,700,800,900,1000])

    msm_ti = get_timescale_data(topologies, top_names, b_values, replicas, coordname, n_timescales)

    raise SystemExit

    if not os.path.exists("plots"):
        os.mkdir("plots")
    os.chdir("plots")
    os.chdir("msm")

    #plot_protein_grid_ti_noavg(msm_ti, topologies, top_names, float_b, replicas, lag_idx, time_idx, True, ti_colors)
    #plot_protein_grid_ti_noavg(msm_ti, topologies, top_names, float_b, replicas, lag_idx, time_idx, False, ti_colors)

    plot_ti_repavg(msm_ti, topologies, top_names, float_b, replicas, lag_idx, time_idx, True, ti_colors)
    #plot_ti_repavg(msm_ti, topologies, top_names, float_b, replicas, lag_idx, time_idx, False, ti_colors)

    #plot_ti_repavg_norm(msm_ti, topologies, top_names, float_b, replicas, lag_idx, time_idx, True, ti_colors)
    #plot_ti_repavg_norm(msm_ti, topologies, top_names, float_b, replicas, lag_idx, time_idx, False, ti_colors)

    os.chdir("..")
    os.chdir("..")

    plt.show()
