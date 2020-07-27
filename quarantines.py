""" File that holds functions/classes for running quarantines.
	Primarily based of off EoN. Will partition this file more later
	as methods get more nuanced

    Throughout: Z := series of I+R (total people who've ever had disease)
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import EoN
import time
import seaborn as sns
import utils as utils

# =====================================================================
# =           OUTPUT DATA STRUCTURE                                   =
# =====================================================================

class TupleSIR:
    """ Class to handle tuples """
    def __init__(self, t, S, I, R, tmax=None):
        self.t = t
        self.S = S
        self.I = I
        self.R = R
        self.Z = np.array([sum(_) for _ in zip(I, R)])
        if tmax is not None:
            self.tmax = tmax 
        else:
            self.tmax = t[-1]

    # =================================================================
    # =           TupleSIR Constructors                               =
    # =================================================================

    @classmethod
    def from_summary(cls, summ, tmax=None):
        return cls(summ.t(), summ.S(), summ.I(), summ.R(), tmax=tmax)

    @classmethod
    def empty(cls):
        return cls(*[np.empty(0) for _ in range(4)], tmax=0.0)


    @classmethod
    def cat(cls, tup_list):
        # Given a list of TupleSIR objects, concat into a single TupleSIR object
        if isinstance(tup_list, cls):
            return tup_list # should just be a single TupleSIR object

        t_lists, s_lists, i_lists, r_lists = [], [], [], []
        running_t = 0
        running_r = 0
        for tup in tup_list:
            # Save in-run data
            t_lists.append(running_t + tup.t)
            s_lists.append(tup.S)
            i_lists.append(tup.I)
            r_lists.append(running_r + tup.R)

            # add end-time-arrays
            t_lists.append(np.array([tup.tmax + running_t]))
            s_lists.append(np.array([tup.S[-1]]))
            i_lists.append(np.array([tup.I[-1]]))
            r_lists.append(np.array([tup.R[-1] + running_r]))

            # update running_r,t
            running_t += tup.tmax
            running_r += tup.I[-1] + tup.R[-1]


        # And concatenate everything in the end
        final_t = np.concatenate(t_lists)
        final_s = np.concatenate(s_lists)
        final_i = np.concatenate(i_lists)
        final_r = np.concatenate(r_lists)

        return TupleSIR(final_t, final_s, final_i, final_r, tmax=running_t)

    # ======  End of TupleSIR Constructors                      =======


    def get_by_str(self, s):
        return getattr(self, s)

    def get_max_I(self):
        return max(self.I)

    def get_final_R(self):
        return self.R[-1]

    def get_peak_I_time(self):
        # Returns (time_of_max_I, max_I)
        return utils.argmax(self.I)

    def plot_single(self, plot_series='IR'):
        plot_sir_counts([self], plot_series=plot_series)

    def plot_I_plus_R(self):
        I_plus_R = [sum(_) for _ in zip(self.I, self.R)]
        fig, ax = plt.subplots(figsize=(8,8))
        ax.plot(self.t, I_plus_R, c='b', alpha=0.5)
        ax.set_xlabel('Time')
        ax.set_ylabel('# I + # R')
        ax.set_title("I+R vs Time")


class AggregateTuple:
    def __init__(self, tup_list):
        """ Placeholder class to hold/access lists of TupleSIR objects """
        self.tup_list = tup_list 

    def __getitem__(self, i):
        return self.tup_list[i]

    def get_max_I(self):
        return utils.mean(self.tup_list, TupleSIR.get_max_I)

    def get_final_R(self):
        return utils.mean(self.tup_list, TupleSIR.get_final_R)

    def get_peak_I_time(self):
        return utils.mean(self.tup_list, TUpleSIR.get_peak_I_time)

    def plot_single(self, plot_series='IR'):
        for serie in plot_series:
            fig, ax = plt.subplots(figsize=(8,8))
            for i, tup in enumerate(self.tup_list):
                ax.plot(tup.t, tup.get_by_str(serie), c='k', alpha=0.5)
            ax.set_xlabel('Time')
            ax.set_ylabel(serie)
            ax.set_title(serie + ' by time')
            ax.legend() 
        


# =========================================================================
# =           BASIC SIR-RUNNING METHODS                                   =
# =========================================================================

def run_until_time(G, tau, gamma, rho, tmax):
    """ Runs basic SIR until time T:
    
    RETURNS:
        G', TupleSIR
        G': new graph object with I,R nodes removed (unless tmax=0.0, then returns identical G)
        TUpleSIR: TupleSIR object with run data
    """
    # handle edge case where tmax == 0:

    if tmax == 0:
        return G, TupleSIR.empty() 
    summary = EoN.fast_SIR(G, tau, gamma, rho=rho, tmax=tmax, return_full_data=True)
    
    summary_dict = utils.invert_dict(summary.get_statuses(time=tmax))
    
    I_nodes = summary_dict.get('I', [])
    R_nodes = summary_dict.get('R', [])

    G = G.copy() 
    G.remove_nodes_from(I_nodes)
    G.remove_nodes_from(R_nodes)
    return G, TupleSIR.from_summary(summary, tmax=tmax)


def run_until_prop_IR(G, tau, gamma, rho, tmax, prop, total_nodes=None):
    """
    Runs SIR model until prop (in [0,1]) fraction of nodes are in I+R states 
    If total_nodes is not None, then the proportion is WRT total_nodes (and not len(G))
    
    RETURNS: G2, status_dict
    G2 is a copy of the graph G with I,R nodes removed 
    status_dict is a dict with keys {t, S, I, R} pointing to the right times 
    """
    if min([tmax, prop]) == 0:
        return G, TupleSIR.empty()

    total_nodes = total_nodes or len(G)
    threshold = total_nodes * prop
    # This has to be slower because we need to run the infection 
    # and then figure out which time to cut things off (ex-post facto)
    summary = EoN.fast_SIR(G, tau, gamma, rho=rho, tmax=tmax, return_full_data=True)

    I, R = summary.I(), summary.R() 
    breakpoint = None
    for i in range(len(I)):
        if I[i] + R[i] >= threshold:
            breakpoint = i 
            break
            
    if breakpoint is not None: # if achieved threshold, modify the graph
        breaktime = summary.t()[breakpoint]
    else:
        breaktime = tmax 
        breakpoint = len(I)
    summary_dict = utils.invert_dict(summary.get_statuses(time=breaktime))
    #print("summary_dict", [(k, len(v)) for k,v in summary_dict.items()])
    I_nodes = summary_dict.get('I', [])
    R_nodes = summary_dict.get('R', [])
    G = G.copy() 
    G.remove_nodes_from(I_nodes)
    G.remove_nodes_from(R_nodes)
    trunc_t = summary.t()[:breakpoint + 1]
    trunc_S = summary.S()[:breakpoint + 1]
    trunc_I = summary.I()[:breakpoint + 1]
    trunc_R = summary.R()[:breakpoint + 1]
    return G, TupleSIR(trunc_t, trunc_S, trunc_I, trunc_R, tmax=breaktime)



# =========================================================================
# =           SINGLE RUNS WITH QUARANTINE(S)                              =
# =========================================================================


def quarantines_by_time(G, tau, gamma, rho, qtimes, tmax, num_iter=1):
    # Runs a model for time tmax, with quarantines specified by qtimes
    # Qtimes is a list of cumulative quarantine times
    tup_list = []

    if not isinstance(qtimes, list):
        qtimes = [qtimes]
    deltas = [qtimes[0]]

    for i in range(len(qtimes) - 1):
        deltas.append(qtimes[i + 1] - qtimes[i])
    deltas.append(tmax - qtimes[-1])
    outputs = []
    original_G = G 
    for iter_num in range(num_iter):
        G = original_G
        tups = []
        for delta in deltas:
            G, tup = run_until_time(G, tau, gamma, rho=rho, tmax=delta)
            tups.append(tup)
        outputs.append(TupleSIR.cat(tups))

    if len(outputs) == 1:
        return outputs[0] 
    else:
        return AggregateTuple(outputs)


def quarantine_by_prop(G, tau, gamma, rho, prop_list, tmax, num_iter=1):
    # SINGLE QUARANTINE ONLY!!!
    # Runs a single quarantine which is initialized when the proportion I,R gets to prop
    # and then runs until tmax afterwards

    if not isinstance(prop_list, list):
        prop_list = [prop_list]

    outputs = [] 
    original_G = G
    for iter_num in range(num_iter):
        G = original_G 
        tups = [] 
        remaining_time = tmax
        for prop in prop_list:
            G, tup = run_until_prop_IR(G, tau, gamma, rho, remaining_time, 
                                       prop, total_nodes=len(original_G))
            remaining_time -= tup.tmax
            tups.append(tup)            
        if remaining_time > 0.0:
            tups.append(run_until_time(G, tau, gamma, rho, remaining_time)[1])
        outputs.append(TupleSIR.cat(tups))

    if len(outputs) == 1:
        return outputs[0] 
    else:
        return AggregateTuple(outputs)


# ======================================================================
# =           PLOTTING METHODS                                         =
# ======================================================================




SERIES_IDX = {'S': 1, 'I': 2, 'R': 3}
def plot_vanilla_run(G, tau, gamma, rho, tmax, series='IR'):
    axlist = []
    for serie in series: 
        fig, ax = plt.subplots(figsize=(8,8))

        select = lambda tup: tup[SERIES_IDX[serie]]
        for i in range(5):
            runtup = EoN.fast_SIR(G, tau, gamma, rho=rho, tmax=tmax)
            ax.plot(runtup[0], select(runtup), c='k', alpha=0.3)
        axlist.append(ax)
    if len(axlist) > 1:
        return axlist
    return axlist[0]


# THINGS TO PLOT 
# X-axis: timestep quarantine was performed at 
# Y-axis1: maximum number of infected at any one time 
# Y-axis2: number of recovered at t=infinity (so number of people who got it )


def plot_single_qs(qrange, G, tau, gamma, rho, maxt):
    # qrange is a list of 
    #step 1, run each simulation:
    sum_ranges = [] 
    for qseries in qrange:
        print("Running sim on series:", qseries)
        sum_ranges.append(quarantine_cycle(G, tau, gamma, rho, qseries, maxt))
        
    
    #step2, triples we care about (quarantine_time, maxI, finalR)
    triples = [(qrange[i], get_max_I(_), get_final_R(_)) for i,_ in enumerate(sum_ranges)]
    
    # step 3, plot both series
    plt.plot(utils.tuple_filter(triples, 0), utils.tuple_filter(triples, 1),
             c='b', alpha=0.5, label='maxI')
    plt.plot(utils.tuple_filter(triples, 0), utils.tuple_filter(triples, 2),
             c='r', alpha=0.5, label='finalR')
    
    return sum_ranges



def plot_sir_counts(tup_list, plot_series='IR', ax=None):
    # Makes a plot for each of the plots in plot_series
    str_dict = {'I': '# Infected', 
                  'R': '# Recovered'}
    for serie in plot_series:
        fig, ax = plt.subplots(figsize=(8,8))
        for i, tup in enumerate(tup_list):
            ax.plot(tup.t, tup.get_by_str(serie), c=utils.c(i), alpha=0.5, label=getattr(tup, 'label', None))
        ax.set_xlabel('Time')
        ax.set_ylabel(str_dict[serie])
        ax.set_title(str_dict[serie] + ' by time')
        ax.legend() 
    


# =====================================================================
# =           2 Quarantines with Grid/Heatmapping                     =
# =====================================================================

def get_quarantine_grid_data_times(G, tau, gamma, rho, tmax, 
                             first_qrange, second_qrange, num_iter=3):
    qpairs = [[q0, q1] for q0 in first_qrange for q1 in second_qrange]
    output_data = {} 
    for q0, q1 in qpairs:
        print("Quarantine @ times t=%s,%s" % (q0, q0 + q1))
        output_data[(q0, q1)] = quarantines_by_time(G, tau, gamma, rho, 
                                                    [q0, q0 + q1], tmax, 
                                                    num_iter=num_iter)
    return output_data 


def get_quarantine_grid_data_props(G, tau, gamma, rho, tmax, first_proprange,
                                   second_proprange, num_iter=3):
    prop_pairs = [[p0, p1] for p0 in first_proprange for p1 in second_proprange]
    output_data = {} 
    for p0, p1 in prop_pairs:
        print("Quarantine @ props p=%s,%s", (p0, p1))
        output_data[(p0, p1)] = quarantine_by_prop(G, tau, gamma, rho, 
                                                   [p0, p1], tmax, 
                                                   num_iter=num_iter)
    return output_data


def process_into_grid(data_dict, func=None):
    if func is None:
        func = lambda x: x 

    # First make a grid from the data:
    keys = sorted(data_dict.keys())
    # make rows:(based on first element of key)
    grid = []
    prev_key = keys[0][0]
    current_row = []
    grid_check = []
    current_row_grid_check = []
    for k in keys:
        if k[0] != prev_key:
            grid.append(current_row)
            current_row = [] 
            prev_key = k[0]
            grid_check.append(current_row_grid_check)
            current_row_grid_check = []
        current_row.append(func(data_dict[k]))
        current_row_grid_check.append(k)
    grid.append(current_row)
    grid_check.append(current_row_grid_check)
    return grid, grid_check



def heatmapify(grid, grid_idxs, title=None):
    """ Builds a heatmap and formats it (except for the title) """
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(grid, ax=ax)
    yticks = [_[0][0] for _ in grid_idxs]
    xticks = [_[1] for _ in grid_idxs[0]]
    ax.set_ylabel("Time of first quarantine")
    ax.set_xlabel("Time of second quarantine (after first)")
    ax.set_yticklabels(yticks)
    ax.set_xticklabels(xticks)
    if title is not None:
        ax.set_title(title)
    return ax


# =========================================================================
# =           Degree map change stuff                                     =
# =========================================================================

""" Second/cleaner attempt at quantifying how degrees change after
    each quarantine 
"""


