import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
cmap = cm.get_cmap("viridis") 

import project_util

def plot_S_Enat_Enon():
    # Plot S(E_non) for each stratum of E_nat.
    plt.figure()
    for i in range(nbins_Enat):
        use_bins = (np.isnan(S_Enat_Enon[i,:]) == False) & (S_Enat_Enon[i,:] > 0)
        if np.sum(use_bins) >= 3:
            S_Enon = S_Enat_Enon[i,use_bins]
            Enon_temp = Enon_mid_bin[use_bins]
            plt.plot(Enon_temp, S_Enon, label="{:.2f}".format(Enat_mid_bin[i]), color=cmap(color_idxs[i]))

    #plt.legend(loc=2)
    plt.xlabel("$E_{non}$ (k$_B$T)")
    plt.ylabel("$S(E_{non})$ (k$_B$)")

    #plt.figure()
    #plt.plot(Enat_mid_bin, Tg_vs_Enat)
    #plt.xlabel("$E_{nat}$")
    #plt.ylabel("$T_g$")
    #plt.title("Glass temperature")
    #plt.show()

    plt.figure()
    plt.pcolormesh(Enat_grid, Enon_grid, S_Enat_Enon_ma, cmap=cmap) # Does this need to be transposed?
    plt.xlabel("$E_{nat}$")
    plt.ylabel("$E_{non}$")
    plt.title("Entropy $S(E_{nat}, E_{non})$")
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    topologies = ["alpha", "beta", "mixed"]
    top_names = [["1imq", "1r69"], ["1fmk", "2akk"], ["1e0g"]]

    topologies = ["alpha"]
    top_names = [["1imq"]]

    #topologies = ["beta"]
    #top_names = [["2akk"]]

    b_values = ["0.01", "0.10", "0.30", "0.50", "0.75", "0.83", "0.92","1.00",
                "1.04", "1.08", "1.12", "1.16", "1.20", "1.25"]

    replicas = range(1,11)

