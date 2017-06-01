import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

import pyemma.msm as msm

import project_util

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
                    if os.path.exists("msm/timescales.npy"):
                        # load in implied timescales
                        lagtimes = np.load("msm/lagtimes.npy")
                        ti = np.load("msm/timescales.npy")
                    else:
                        ti = np.nan*np.zeros((16,4))

                    msm_ti_rep.append(ti)
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

if __name__ == "__main__":
    colorcycle = ["#5DA5DA", "#FAA43A", "#60BD68", "#F17CB0", "#B2912F", "#B276B2", "#DECF3F", "#F15854", "#4D4D4D"]

    #topologies = ["alpha", "beta"]
    #top_names = [["1imq"], ["2akk"]]

    #topologies = ["alpha"]
    #top_names = [["1imq"]]

    topologies = ["alpha", "beta", "mixed"]
    symbols = ["^", "o", "s"]
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    color_codes = [["#5DA5DA", "#F17CB0"], ["#60BD68", "#B2912F"], ["#FAA43A"]]

    #topologies = ["alpha", "beta"]
    #top_names = [["1r69", "1imq"], ["1fmk"]]

    topologies = ["beta"]
    top_names = [["1fmk"]]

    #b_values = ["0.00", "0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
    #            "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]
    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]

    #replicas = [1,2,3,4,5,6,7,8,9]
    replicas = range(1, 11)

    coordname = "Qtanh_0_05"
    datapath = "msm/Neff.dat"
    semilogy = True
    n_timescales = 4
    #lagtime = 200
    #lag_idx = 7
    ti_colors = ['r','g','b','salmon', 'darkgreen']

    savefig = False

    lags = np.array([1,2,5,10,20,50,100,200,300,400,500,600,700,800,900,1000])
    lagtime_idx = 7 # 200

    msm_ti = get_timescale_data(topologies, top_names, b_values, replicas, coordname, n_timescales)

    Neff_data = project_util.Dataset(topologies, top_names, b_values, replicas, datapath)

    if not os.path.exists("plots"):
        os.mkdir("plots")
    os.chdir("plots")
    os.chdir("msm")

    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            fig, axes = plt.subplots(len(float_b), len(replicas),
                figsize=(20,20), sharex=True, sharey='row')

            for b in range(len(float_b)):
                for rep in range(len(replicas)):
                    ax = axes[b, rep]


                    # plot implied timescales for 
                    rep_its = msm_ti[t][n][b][rep]
                    if not np.all(np.isnan(rep_its)):
                        for i in range(len(rep_its[0])):
                            # check for NaNs
                            if not np.any(np.isnan(rep_its[:,i])):
                                if i == 0:
                                    ax.annotate("$N_{{eff}} = {}$".format(Neff_data.data[t][n][b][rep]), 
                                            fontsize=8, xy=(0,0), xytext=(0.45, 0.1),
                                            textcoords="axes fraction", xycoords="axes fraction")
                                ax.plot(lags, rep_its[:,i])
                                ax.fill_between(lags, np.ones(len(lags)), lags, facecolor='gray', alpha=1)
                                ax.plot(lags, lags, 'k', lw=1)
                                if semilogy:
                                    ax.semilogy()

                    
                        if b == (len(float_b) - 1):
                            ax.set_xticklabels([0,200,400,600,800], rotation="vertical", fontsize=8)

                        if rep == 0:
                            for tick in ax.yaxis.get_major_ticks():
                                tick.label.set_fontsize(8) 
                    else:
                        for tick in ax.yaxis.get_major_ticks():
                            tick.label.set_visible(False) 
                        for tick in ax.xaxis.get_major_ticks():
                            tick.label.set_visible(False) 

            plt.subplots_adjust(wspace=0)

            axwidth = 1./float(len(replicas))
            axheight = 1./float(len(float_b))

            big_ax = fig.add_subplot(111)
            big_ax.grid(False)
            big_ax.set_axis_bgcolor('none')
            big_ax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
            big_ax.set_title(top_names[t][n], fontsize=20)
            big_ax.annotate("b-values", fontsize=18, rotation="vertical",
                xy=(0,0), xytext=(-axwidth, 0.45),
                textcoords="axes fraction", xycoords="axes fraction")

            big_ax.annotate("Replicas", fontsize=18,
                xy=(0,0), xytext=(0.45, -0.65*axheight),
                textcoords="axes fraction", xycoords="axes fraction")

            for x in range(len(replicas)):
                big_ax.annotate(str(x + 1), xy=(0,0), fontsize=16,
                    xytext=(0.5*axwidth + x*axwidth, -0.5*axheight),
                    textcoords="axes fraction", xycoords="axes fraction")

            for x in range(len(float_b)):
                big_ax.annotate("{:.2f}".format(float_b[-(x + 1)]), xy=(0,0),
                    fontsize=14, 
                    xytext=(-0.75*axwidth, 0.5*axheight + x*axheight),
                    textcoords="axes fraction", xycoords="axes fraction")
            if savefig:
                if semilogy:
                    fig.savefig("{}_its_vs_lag_semilog.pdf".format(top_names[t][n]))
                else:
                    fig.savefig("{}_its_vs_lag.pdf".format(top_names[t][n]))
    os.chdir("..")
    os.chdir("..")
    plt.show()
