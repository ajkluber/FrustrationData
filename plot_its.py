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

                    if os.path.exists("msm/msm.pkl"):
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

    #topologies = ["beta"]
    #top_names = [["2akk"]]

    #b_values = ["0.00", "0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
    #            "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]
    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]

    replicas = range(1,11)

    coordname = "Qtanh_0_05"
    n_timescales = 10
    lagtime = 200
    lag_idx = 7
    ti_colors = ['r','g','b','salmon', 'darkgreen']
    lags = np.array([1,2,5,10,20,50,100,200,300,400,500,600,700,800,900,1000])

    msm_ti = get_timescale_data(topologies, top_names, b_values, replicas, coordname, n_timescales)

    raise SystemExit

    if not os.path.exists("plots"):
        os.mkdir("plots")
    os.chdir("plots")
    os.chdir("msm")

#    #plot_timescales_all_replicas(msm_ti, topologies, top_names): 
#   # ToDo: Update to account for lagtime dimension
#    for t in range(len(topologies)):
#        for n in range(len(top_names[t])):
#            plt.figure()
#            plt.title(top_names[t][n])
#            for i in range(n_timescales): 
#                # plot all replicas
#                first_time = True
#                b_rep_ti = [ x[:,i] for x in msm_ti[t][n] if len(x) > 0 ]
#                b_rep_idxs = [ x for x in range(len(msm_ti[t][n])) if len(msm_ti[t][n][x]) > 0 ]
#                for b in range(len(b_rep_ti)):
#                    for r in range(len(b_rep_ti[b])):
#                        if first_time:
#                            plt.plot(float_b[b_rep_idxs[b]], b_rep_ti[b][r], marker='o', label="$t_"+str(i)+"$", color=ti_colors[i])
#                            first_time = False
#                        else:
#                            plt.plot(float_b[b_rep_idxs[b]], b_rep_ti[b][r], marker='o', color=ti_colors[i])
#        
#            plt.legend(loc=2)
#            plt.show()
#

    # plot normalized disorder-average timescales.
#    for i in range(n_timescales): 
#        plt.figure()
#        for t in range(len(topologies)):
#            for n in range(len(top_names[t])):
#                first_time = True
#                for b in range(len(b_value_has_data)):
#                    # Extract lagtime at this b value
#                    n_reps = len(msm_ti[t][n][b])
#                    name_ti = np.array([ msm_ti[t][n][b][rep][lag_idx,i] for rep in range(n_reps) ])
#                    avg_ti_unfrust = np.mean(np.array([ msm_ti[t][n][0][rep][lag_idx,i] for rep in range(n_reps) ]))
#
#                    #b_value_has_data = np.array([ (len(b_rep_ti) > 0) for x in msm_ti[t][n] ])
#                    #if b_value_has_data[b]:
#                    if len(name_ti) > 0:
#                        #avg_ti_unfrust = np.mean(msm_ti[t][n][0][0][:,i])
#                        # plot disorder average
#                        if first_time:
#                            plt.plot(float_b[b], np.mean(name_ti)/avg_ti_unfrust, marker=symbols[t], color=color_codes[t][n], label=top_names[t][n])
#                            first_time = False
#                        else:
#                            plt.plot(float_b[b], np.mean(name_ti)/avg_ti_unfrust, marker=symbols[t], color=color_codes[t][n])
#        
#        plt.legend(loc=2)
#        ymin, ymax = plt.ylim()
#        plt.ylim(0, ymax)
#        plt.title("MSM $t_" + str(i + 1) + "$")
#        plt.ylabel("Timescales $t_" + str(i + 1) + "$ ($t_" + str(i + 1) + "^{b=0}$)")
#        plt.xlabel("Frustration $b$")
#        plt.savefig("repavg_msm_t_" + str(i + 1) + ".png")
#        plt.savefig("repavg_msm_t_" + str(i + 1) + ".pdf")
#    plt.show()


    # plot implied timescales as a function of lagtime for all proteins and
    # replicas. Each b-value is one figure. 
    for b in range(len(float_b)):
        plt.figure()
        for t in range(len(topologies)):
            for n in range(len(top_names[t])):
                first_time = True
                # plot all replica implied timescales 
                n_reps = len(msm_ti[t][n][b])
                for rep in range(n_reps):
                    rep_its = msm_ti[t][n][b][rep]
                    for i in range(len(rep_its[0])):
                        if first_time:
                            plt.plot(lags, rep_its[:,i], color=color_codes[t][n], label=top_names[t][n])
                            first_time = False
                        else:
                            plt.plot(lags, rep_its[:,i], color=color_codes[t][n])

        plt.fill_between(lags, np.ones(len(lags)), lags, facecolor='gray', alpha=1)
        plt.plot(lags, lags, 'k', lw=2)
        plt.legend(loc=4)
        plt.title("Frustration $b = {:.2f}$".format(float_b[b]))
        plt.ylabel("Implied Timescales")
        plt.xlabel("Lagtime $\\tau$")
        plt.semilogy()
        #plt.savefig("allreps_msm_its_b_{:.2f}.png".format(float_b[b]))
        #plt.savefig("allreps_msm_its_b_{:.2f}.png".format(float_b[b]))
        #break
        plt.show()


#    # plot disorder-average timescales.
#    for t in range(len(topologies)):
#        for n in range(len(top_names[t])):
#            #fig, axes = plt.subplots(2,2,sharex=True, figsize=(12,12))
#            #fig.suptitle(top_names[t][n])
#            plt.figure()
#            for i in range(n_timescales): 
#                #ax = axes[i / 2, i % 2]
#                #ax.annotate("$t_" + str(i) + "$", xy=(0,0), xytext=(0.2,0.5),
#                #        textcoords="axes fraction", xycoords="axes fraction", fontsize=20)
#
#                b_value_has_data = np.array([ (len(x) > 0) for x in msm_ti[t][n] ])
#                first_time = True
#                
#                for b in range(len(b_value_has_data)):
#                    if b_value_has_data[b]:
#                        avg_ti_unfrust = np.mean(msm_ti[t][n][0][:,i])
#
#                        # plot disorder average
#                        name_ti = msm_ti[t][n][b][:,i] 
#                        if first_time:
#                            #ax.plot(float_b[b], np.mean(name_ti), marker='o', color=ti_colors[i], label="$t_{}$".format(i))
#                            plt.plot(float_b[b], np.mean(name_ti), marker=symbols[t], color=color_codes[t][n])
#                            first_time = False
#                        else:
#                            #ax.plot(float_b[b], np.mean(name_ti), marker='o', color=ti_colors[i])
#                            plt.plot(float_b[b], np.mean(name_ti), marker=symbols[t], color=color_codes[t][n])
#            #plt.legend(loc=2, fontsize=16)
#            plt.title(top_names[t][n] + " " + topologies[t] + " N = " + str(prot_sizes[t][n]))
#            plt.xlabel("Frustration $b$")
#            plt.ylabel("MSM timescale $t_i$")
#            plt.savefig(top_names[t][n] + "_msm_its.png")
#            plt.savefig(top_names[t][n] + "_msm_its.pdf")
#    plt.show()

    os.chdir("..")
    os.chdir("..")
