import os
import glob
import pickle
import numpy as np

from scipy.interpolate import interp1d

##########################################################
# CLASSES TO HELP ORGANIZE THE DATASET
##########################################################
class Dataset(object):
    """Gather data that is one-value for each replica"""

    def __init__(self, topologies, top_names, b_values, replicas, datapath,
            skip=[], T_used=False, empty=False, dont_avg=False):

        self.topologies = topologies 
        self.top_names = top_names
        self.b_values = b_values
        self.replicas = replicas
        self.datapath = datapath
        self.skip = skip
        self.T_used = T_used
        self.empty = empty
        self.dont_avg = dont_avg
        self.has_avg = True

        self._get_data()

    def get_replica_data(self):

        if self.T_used:
            dir = os.path.dirname(self.datapath)
            with open(dir + "/T_used.dat", "r") as fin:
                T_used = float(fin.read())
            with open(self.datapath.format(T_used), "r") as fin:
                data = float(fin.read())
        else:
            with open(self.datapath, "r") as fin:
                data = float(fin.read())
        return data

    def _get_data(self):

        all_data = []
        for t in range(len(self.topologies)):
            os.chdir(self.topologies[t])
            names = self.top_names[t]
            data_top = []
            for name in names:
                data_name = []
                if not self.empty:
                    print name
                os.chdir(name)
                for b in self.b_values:
                    if not self.empty:
                        print "getting b-value: ", b
                    data_rep = []
                    for rep in self.replicas:
                        if not os.path.exists("b_{}/replica_{}".format(b, rep)):
                            data_rep.append(np.nan)
                        else:
                            os.chdir("b_{}/replica_{}".format(b, rep))
                            cwd = os.getcwd()
                            if cwd in self.skip: 
                                os.chdir("../..")
                                continue
                            if not self.empty:
                                if os.path.exists(self.datapath):
                                    data = self.get_replica_data()
                                    data_rep.append(data)
                                else:
                                    data_rep.append(np.nan)
                            os.chdir("../..")
                            if b == "0.00":
                                # if unfrustrated go-model only run one replica.
                                break
                    if self.empty:
                        data_name.append(np.nan*np.ones(len(self.replicas)))
                    else:
                        data_name.append(np.array(data_rep))
                os.chdir("..")
                data_top.append(data_name)
            all_data.append(data_top)
            os.chdir("..")
        self.data = all_data
        if not self.dont_avg:
            self._calc_repavg()

    def _calc_repavg(self):
        """ calculate the replicas average data"""

        all_data_avg = []
        all_data_std = []
        for t in range(len(self.topologies)):
            data_top_avg = []
            data_top_std = []
            for n in range(len(self.top_names[t])):
                data_name_avg = []
                data_name_std = []
                for b in range(len(self.b_values)):
                    y = np.array(self.data[t][n][b])
                    good = ~np.isnan(y)
                    if np.sum(good) > 0:
                        data_name_avg.append(np.mean(y[good]))
                        data_name_std.append(np.std(y[good])/np.sqrt(float(np.sum(good))))
                    else:
                        data_name_avg.append(np.nan)
                        data_name_std.append(np.nan)

                data_top_avg.append(np.array(data_name_avg))
                data_top_std.append(np.array(data_name_std))

            all_data_avg.append(data_top_avg)
            all_data_std.append(data_top_std)

        self.avgdata = all_data_avg
        self.stddata = all_data_std
        self.has_avg = True

class DatasetXvsY(object):
    """Gather data that is Y vs X for each replica"""

    def __init__(self, topologies, top_names, b_values, replicas, datapath_x,
            datapath_y, skip=[], T_used=False, empty=False, calc_repavg=True):

        self.topologies = topologies 
        self.top_names = top_names
        self.b_values = b_values
        self.replicas = replicas
        self.datapath_x = datapath_x
        self.datapath_y = datapath_y
        self.skip = skip
        self.T_used = T_used
        self.empty = empty

        self._get_xydata()
        if calc_repavg:
            self._calc_repavg()

    def get_replica_data(self):
        if self.datapath_x.endswith(".npy"):
            load_func_x = np.load
        else:
            load_func_x = np.loadtxt

        if self.datapath_y.endswith(".npy"):
            load_func_y = np.load
        else:
            load_func_y = np.loadtxt

        if self.T_used:
            dir = os.path.dirname(self.datapath_x)
            with open(dir + "/T_used.dat", "r") as fin:
                T_used = float(fin.read())
            data_x = load_func_x(self.datapath_x.format(T_used))
            data_y = load_func_y(self.datapath_y.format(T_used))
        else:
            data_x = load_func_x(self.datapath_x)
            data_y = load_func_y(self.datapath_y)
        return data_x, data_y

    def _get_xydata(self):

        all_data_y = []
        all_data_x = []
        for t in range(len(self.topologies)):
            os.chdir(self.topologies[t])
            names = self.top_names[t]
            data_x_top = []
            data_y_top = []
            for name in names:
                if not self.empty:
                    print name
                data_x_name = []
                data_y_name = []
                os.chdir(name)
                self.unfolded_Q = None
                self.n_res = None
                self.non_pairs = None
                for b in self.b_values:
                    if not self.empty:
                        print "getting b-value: ", b
                    data_x_rep = []
                    data_y_rep = []
                    if not self.empty:
                        for rep in self.replicas:
                            if os.path.exists("b_{}/replica_{}".format(b, rep)):
                                os.chdir("b_{}/replica_{}".format(b, rep))
                                cwd = os.getcwd()
                                if cwd in self.skip: 
                                    os.chdir("../..")
                                    continue
                                T_used = get_T_used()
                                if self.T_used and T_used:
                                    with open("Qtanh_0_05_profile/T_used.dat", "r") as fin:
                                        T_used = float(fin.read())
                                    if os.path.exists(self.datapath_x.format(T_used)) and \
                                            os.path.exists(self.datapath_y.format(T_used)):
                                        data_x, data_y = self.get_replica_data()
                                        data_x_rep.append(data_x)
                                        data_y_rep.append(data_y)
                                else:
                                    if os.path.exists(self.datapath_x) and \
                                            os.path.exists(self.datapath_y):
                                        data_x, data_y = self.get_replica_data()
                                        data_x_rep.append(data_x)
                                        data_y_rep.append(data_y)
                                os.chdir("../..")

                            if b == "0.00":
                                # if unfrustrated go-model only run one replica.
                                break
                    data_x_name.append(data_x_rep)
                    data_y_name.append(data_y_rep)
                os.chdir("..")
                data_x_top.append(data_x_name)
                data_y_top.append(data_y_name)
            all_data_x.append(data_x_top)
            all_data_y.append(data_y_top)
            os.chdir("..")
        self.xdata = all_data_x
        self.ydata = all_data_y

    def _calc_repavg(self):
        """ calculate the replicas average data"""

        all_data_y_avg = []
        all_data_x_avg = []
        for t in range(len(self.topologies)):
            data_x_top_avg = []
            data_y_top_avg = []
            for n in range(len(self.top_names[t])):
                data_x_name_avg = []
                data_y_name_avg = []
                for b in range(len(self.b_values)):
                    data_x_b_avg, data_y_b_avg = self._get_xy_avg(self.xdata[t][n][b], self.ydata[t][n][b])
                    data_x_name_avg.append(data_x_b_avg)
                    data_y_name_avg.append(data_y_b_avg)

                data_x_top_avg.append(data_x_name_avg)
                data_y_top_avg.append(data_y_name_avg)

            all_data_x_avg.append(data_x_top_avg)
            all_data_y_avg.append(data_y_top_avg)

        self.avgxdata = all_data_x_avg
        self.avgydata = all_data_y_avg

    def _get_xy_avg(self, rep_xdata, rep_ydata):
        """Take the replica average"""

        # calculate disordered-averaged profile
        data_x_name_avg, data_y_name_avg = [], [] 
        good = np.array([ np.all(~np.isnan(rep_ydata[r])) for r in range(len(rep_ydata)) ])
        if (len(rep_xdata) > 0) and (np.sum(good) > 0):
            # only average over replicas
            if hasattr(rep_ydata[0], "mask"):
                # determine a shared domain of data
                xmin = np.max([ np.min(rep_xdata[r][~ rep_ydata[r].mask]) for r in range(len(rep_xdata)) if good[r]])
                xmax = np.min([ np.max(rep_xdata[r][~ rep_ydata[r].mask]) for r in range(len(rep_xdata)) if good[r]])
                xinterp = np.linspace(xmin, xmax, 500)
                # only interpolate with non-masked values
                yinterp = np.array([ interp1d(rep_xdata[r][~ rep_ydata[r].mask], rep_ydata[r][~ rep_ydata[r].mask], kind="cubic")(xinterp) for r in range(len(rep_xdata)) if good[r]])
            else:
                xmin = np.max([ np.min(rep_xdata[r]) for r in range(len(rep_xdata)) if good[r]])
                xmax = np.min([ np.max(rep_xdata[r]) for r in range(len(rep_xdata)) if good[r]])
                xinterp = np.linspace(xmin, xmax, 500)
                yinterp = np.array([ interp1d(rep_xdata[r], rep_ydata[r], kind="cubic")(xinterp) for r in range(len(rep_xdata)) if good[r]])
            data_x_name_avg = xinterp
            data_y_name_avg = np.mean(yinterp, axis=0) 
        return data_x_name_avg, data_y_name_avg

##########################################################
# LITTLE HELPER FUNCTIONS
##########################################################
def save_avg_plot_data(dataset, dataset_saveas):
    """Save the parameter-averaged data in central location"""

    if not dataset.has_avg:
        dataset._calc_repavg()

    yavg = dataset.avgdata
    yerr = dataset.stddata
    bvals = np.array([ float(x) for x in dataset.b_values ])

    if not os.path.exists("plots/plot_data"):
        os.mkdir("plots/plot_data")

    data_pkl = {"names":dataset.top_names, "bvals":bvals, "avg":yavg, "err":yerr}

    with open("plots/plot_data/{}.pkl".format(dataset_saveas), "wb") as fout:
        pickle.dump(data_pkl, fout, pickle.HIGHEST_PROTOCOL)


def useable_data(x):
    """ """ 
    good = np.zeros(len(x), bool)
    for i in range(len(x)):
        if x[i] is None:
            good[i] = False
        elif x[i] == []:
            good[i] = False
        elif np.isnan(x[i]):
            good[i] = False
        else:
            good[i] = True
    return good


def get_n_native_pairs(name):
    if os.path.exists(name + ".contacts"):
        n_native_pairs = len(np.loadtxt(name + ".contacts"))
    elif os.path.exists(name + ".ini"):
        with open(name + ".ini", "r") as fin:
            n_native_pairs = int([ x for x in fin.readlines() if x.startswith("n_native_pairs") ][0].split()[-1])
    else:
        return False
    return n_native_pairs

def get_T_used():
    if os.path.exists("Qtanh_0_05_profile/T_used.dat"):
        with open("Qtanh_0_05_profile/T_used.dat", "r") as fin:
            T_used = float(fin.read())
    else:
        T_used = False

    return T_used

def ISdir_exist():
    T = get_T_used()
    exist = np.all(np.array([ os.path.exists(
        'T_{:.2f}_{}/inherent_structures'.format(T, x)) for x in [1,2,3] ]))
    return exist

def engIScat_exist():
    T = get_T_used()
    cwd = os.getcwd()
    exist = []
    for j in range(1,4):
        os.chdir("T_{:.2f}_{}/inherent_structures".format(T, j))
        exist.append(os.path.exists('Ebackbone.npy'))
        exist.append(os.path.exists('Enat.npy'))
        exist.append(os.path.exists('Enon.npy'))
        exist.append(os.path.exists('Ebackbone_thm.npy'))
        exist.append(os.path.exists('Enat_thm.npy'))
        exist.append(os.path.exists('Enon_thm.npy'))
        os.chdir(cwd)

    return np.all(exist) 

def engISrank_exist():
    T = get_T_used()
    cwd = os.getcwd()
    exist = []
    for j in range(1,4):
        os.chdir("T_{:.2f}_{}/inherent_structures".format(T, j))
        size = len(glob.glob("rank_*"))
        exist.append(np.all([ os.path.exists("rank_{}/Ebackbone.npy".format(x)) for x in range(size)]))
        exist.append(np.all([ os.path.exists("rank_{}/Enat.npy".format(x)) for x in range(size)]))
        exist.append(np.all([ os.path.exists("rank_{}/Enon.npy".format(x)) for x in range(size)]))
        exist.append(np.all([ os.path.exists("rank_{}/Ebackbone_thm.npy".format(x)) for x in range(size)]))
        exist.append(np.all([ os.path.exists("rank_{}/Enat_thm.npy".format(x)) for x in range(size)]))
        exist.append(np.all([ os.path.exists("rank_{}/Enon_thm.npy".format(x)) for x in range(size)]))
        os.chdir(cwd)

    return np.all(exist) 

def qtanh_exist():
    T = get_T_used()
    qtanh_files_exist = np.all(np.array([ os.path.exists(
        'T_{:.2f}_{}/Qtanh_0_05.npy'.format(T, x)) for x in [1,2,3] ]))
    return qtanh_files_exist

def trajIScat_exist():
    T = get_T_used()
    cwd = os.getcwd()
    exist = np.all([ os.path.exists("T_{:.2f}_{}/inherent_structures/traj.xtc".format(T, j)) for j in [1, 2, 3]])
    return exist

def trajISrank_exist():
    T = get_T_used()
    cwd = os.getcwd()
    exist = []
    for j in range(1,4):
        os.chdir("T_{:.2f}_{}/inherent_structures".format(T, j))
        size = len(glob.glob("rank_*"))
        exist.append(np.all([ os.path.exists("rank_{}/all_frames.xtc".format(x)) for x in range(size)]))
        os.chdir(cwd)

    return np.all(exist) 

def trajTf_exist():
    T = get_T_used()
    traj_files_exist = np.all(np.array([ os.path.exists(
        'T_{:.2f}_{}/traj.xtc'.format(T, x)) for x in [1,2,3] ]))
    return traj_files_exist

