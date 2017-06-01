import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

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
    float_b = np.array([ float(x) for x in b_values ])

    replicas = range(1,11)

    plotstyle = project_plotter.global_plot_style()

    Tf_plotspecs = {"xlabel":"Frustration ($b$)", "ylabel":"Temperature",
            "legend_loc":6, "saveas":None}
    Tf_plotspecs.update(plotstyle)

    Tf_data = project_util.Dataset(topologies, top_names, b_values, replicas, "Qtanh_0_05_profile/T_used.dat")
    Tf_data._calc_repavg()

    model2 = lambda x, a: np.sqrt(1 - a*(x**2))
    model4 = lambda x, a: np.sqrt(1 - a*(x**4))

    coeffs2 = []
    coeffs4 = []
    for t in range(len(topologies)):
        for n in range(len(top_names[t])):
            Tf0 = Tf_data.avgdata[t][n][0]
            y = Tf_data.avgdata[t][n]
            dy = Tf_data.stddata[t][n]
            good = ~np.isnan(y)
            y = y[good]
            dy = dy[good]
            b = float_b[good]

            popt2, pcov = scipy.optimize.curve_fit(model2, b, y/Tf0, p0=(0.0001))

            popt4, pcov = scipy.optimize.curve_fit(model4, b, y/Tf0, p0=(0.0001))
            #print popt2, popt4
            coeffs2.append(popt2[0])
            coeffs4.append(popt4[0])
            
            plt.errorbar(b, y, yerr=1.96*dy, marker=plotstyle["markers"][t], ls="None",
                    color=plotstyle["color"][t][n])

            x = np.linspace(0, max(b), 100)

            if (t == 0) and (n == 0):
                plt.plot(x, Tf0*model4(x, popt4[0]), ls='-',
                        color="k", label=r"$T_f^{b=0}\sqrt{1 - c_4b^4}$")

                plt.plot(x, Tf0*model2(x, popt2[0]), ls='--', 
                        color="k", label=r"$T_f^{b=0}\sqrt{1 - c_2b^2}$")

            plt.plot(x, Tf0*model4(x, popt4[0]), ls='-',
                    color=plotstyle["color"][t][n])

            plt.plot(x, Tf0*model2(x, popt2[0]), ls='--', 
                    color=plotstyle["color"][t][n])

    np.savetxt("plots/Tf_fit_2_coeffs.dat", np.array(coeffs2))
    np.savetxt("plots/Tf_fit_4_coeffs.dat", np.array(coeffs4))

    plt.legend(loc=3)
    plt.xlabel("Frustration ($b$)")
    plt.ylabel(r"$T_f$ (K)")
    plt.savefig("plots/Tf_theory_fit_b2_vs_b4.pdf")
    plt.savefig("plots/Tf_theory_fit_b2_vs_b4.png")
    plt.show()
