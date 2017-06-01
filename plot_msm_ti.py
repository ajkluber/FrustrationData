import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

import pyemma.msm as msm

def get_timescale_data(topologies, top_names, b_values, float_b, replicas, coordname, n_timescales, lagtime):
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
                    if os.path.exists("msm/msm_ti_{}.npy".format(lagtime)):
                        msm_ti_rep.append(np.load("msm/msm_ti_{}.npy".format(lagtime)))
                    else:
                        if os.path.exists("msm/msm.pkl"):
                            try:
                                with open("msm/msm.pkl", "rb") as fhandle:
                                    msm_info = pickle.load(fhandle)
                                    model = msm.markov_model(msm_info[str(lagtime)])
                                    print lagtime*model.timescales()[:n_timescales]
                                    temp_timescales = lagtime*model.timescales()[:n_timescales]
                                    msm_ti_rep.append(temp_timescales)
                                    np.save("msm/msm_ti_{}.npy".format(lagtime), temp_timescales)
                            except:
                                print "unsuccessful:", os.getcwd()
                    os.chdir("../..")
                    if b == "0.00":
                        # if unfrustrated go-model only run one replica.
                        break
                msm_ti_name.append(np.array(msm_ti_rep))
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
    top_names = [["1r69", "1imq"], ["1fmk", "2akk"], ["1e0g"]]
    prot_sizes = [[63, 86], [58, 74], [48]]
    symbols = ["^", "o", "s"]
    color_codes = [["#5DA5DA", "#F17CB0"], ["#60BD68", "#B2912F"], ["#FAA43A"]]

    #topologies = ["beta"]
    #top_names = [["2akk"]]

    #b_values = ["0.00", "0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
    #            "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]
    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92",
                "1.00", "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    float_b = [ float(x) for x in b_values ]

    replicas = [1,2,3,4,5,6,7,8,9]

    coordname = "Qtanh_0_05"
    n_timescales = 4
    lagtime = 200
    ti_colors = ['r','g','b','salmon', 'darkgreen']

    msm_ti = get_timescale_data(topologies, top_names, b_values, float_b, replicas, coordname, n_timescales, lagtime)


    if not os.path.exists("plots"):
        os.mkdir("plots")
    os.chdir("plots")
    os.chdir("msm")

#    #plot_timescales_all_replicas(msm_ti, topologies, top_names):
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
#    plt.figure()
#    for t in range(len(topologies)):
#        for n in range(len(top_names[t])):
#
#            first_time = True
#            for i in range(n_timescales): 
#                b_value_has_data = np.array([ (len(x) > 0) for x in msm_ti[t][n] ])
#                
#                for b in range(len(b_value_has_data)):
#                    if b_value_has_data[b]:
#                        avg_ti_unfrust = np.mean(msm_ti[t][n][0][:,i])
#
#                        # plot disorder average
#                        name_ti = msm_ti[t][n][b][:,i] 
#                        if first_time:
#                            plt.plot(float_b[b], np.mean(name_ti)/avg_ti_unfrust, ms=5, marker=symbols[t], color=ti_colors[i], label=top_names[t][n])
#                            first_time = False
#                        else:
#                            plt.plot(float_b[b], np.mean(name_ti)/avg_ti_unfrust, ms=5, marker=symbols[t], color=ti_colors[i])
#        
#    plt.legend(loc=2)
#    plt.title(" ")
#    plt.ylabel("$t_i$ / $t_i^{b=0}$")
#    plt.xlabel("Frustration $b$")
#    plt.savefig("all_msm_t_i.png")
#    plt.savefig("all_msm_t_i.pdf")
#    plt.show()

    # plot normalized disorder-average timescales.
    for i in range(n_timescales): 
        plt.figure()
        for t in range(len(topologies)):
            for n in range(len(top_names[t])):

                b_value_has_data = np.array([ (len(x) > 0) for x in msm_ti[t][n] ])
                first_time = True
                
                for b in range(len(b_value_has_data)):
                    if b_value_has_data[b]:
                        avg_ti_unfrust = np.mean(msm_ti[t][n][0][:,i])

                        # plot disorder average
                        name_ti = msm_ti[t][n][b][:,i] 
                        if first_time:
                            plt.plot(float_b[b], np.mean(name_ti)/avg_ti_unfrust, ms=8, marker=symbols[t], color=color_codes[t][n], label=top_names[t][n])
                            first_time = False
                        else:
                            plt.plot(float_b[b], np.mean(name_ti)/avg_ti_unfrust, ms=8, marker=symbols[t], color=color_codes[t][n])
        
        plt.legend(loc=2)
        plt.title(" ")
        if i == 0:
            plt.ylabel("$t_{fold}$ / $t_{fold}^{b=0}$")
        else:
            plt.ylabel("$t_{0}$ / $t_{0}^{{b=0}}$".format(i + 1))
        plt.xlabel("Frustration $b$")
        plt.savefig("all_msm_t_" + str(i + 1) + ".png")
        plt.savefig("all_msm_t_" + str(i + 1) + ".pdf")
    plt.show()


    # plot disorder-average timescales.
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            #fig, axes = plt.subplots(2,2,sharex=True, figsize=(12,12))
            #fig.suptitle(top_names[t][n])
            plt.figure()
            for i in range(n_timescales): 
                #ax = axes[i / 2, i % 2]
                #ax.annotate("$t_" + str(i) + "$", xy=(0,0), xytext=(0.2,0.5),
                #        textcoords="axes fraction", xycoords="axes fraction", fontsize=20)

                b_value_has_data = np.array([ (len(x) > 0) for x in msm_ti[t][n] ])
                first_time = True
                
                for b in range(len(b_value_has_data)):
                    if b_value_has_data[b]:
                        avg_ti_unfrust = np.mean(msm_ti[t][n][0][:,i])

                        # plot disorder average
                        name_ti = msm_ti[t][n][b][:,i] 
                        if first_time:
                            #ax.plot(float_b[b], np.mean(name_ti), marker='o', color=ti_colors[i], label="$t_{}$".format(i))
                            plt.plot(float_b[b], np.mean(name_ti), ms=8, marker=symbols[t], color=color_codes[t][n])
                            first_time = False
                        else:
                            #ax.plot(float_b[b], np.mean(name_ti), marker='o', color=ti_colors[i])
                            plt.plot(float_b[b], np.mean(name_ti), ms=8, marker=symbols[t], color=color_codes[t][n])
            #plt.legend(loc=2, fontsize=16)
            plt.title(top_names[t][n] + " " + topologies[t] + " N = " + str(prot_sizes[t][n]))
            plt.xlabel("Frustration $b$")
            plt.ylabel("MSM timescale $t_i$")
            plt.savefig(top_names[t][n] + "_msm_its.png")
            plt.savefig(top_names[t][n] + "_msm_its.pdf")
    plt.show()

    os.chdir("..")
    os.chdir("..")
