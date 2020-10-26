""" File that holds functions/classes for running quarantines.
	Primarily based of off EoN. Will partition this file more later
	as methods get more nuanced

    Throughout: Z := series of I+R (total people who've ever had disease)
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import EoNlocal as EoN # Locally modified version of EoN
import time
import seaborn as sns
import utils as utils
import pickle 
from bson.binary import Binary

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
            if tup.tmax == 0:
                continue

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

    @classmethod
    def from_dict(cls, data_dict):
        tmax = data_dict.get('tmax', data_dict['t'][-1])
        return cls(data_dict['t'], data_dict['S'], data_dict['I'], 
                   data_dict['R'], tmax=tmax)

    def to_dict(self):
        return {'t': self.t, 
                'S': self.S, 
                'I': self.I, 
                'R': self.R, 
                'tmax': self.tmax}
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

    def get_peak_width(self):
        """ Gets the "width" of the peak (assuming this has a single peak)
            by finding the peak I and considering the time-gap between 
            the times the peak hits half-peak
        """
        peak_i_idx, peak_i = max(enumerate(self.I), key=lambda p: p[1])
        halfpeak = peak_i / 2
        # Get first time hitting half:
        pre_idx = peak_i_idx
        post_idx = peak_i_idx
        while self.I[pre_idx] > halfpeak:
            pre_idx -= 1

        while self.I[post_idx] > halfpeak:
            post_idx += 1

        return self.t[post_idx] - self.t[pre_idx]


    def get_all_peak_widths(self):
        """ Gets a list of peak widths (in terms of time) """

        # First find all the locally maximal indices 
        # Very structured signal: 
        # Can identify the indices when quarantine was applied 
        quarantine_idxs = [] 
        for i in range(len(self.I) - 1):
            if self.I[i + 1] < self.I[i] - 1:
                quarantine_idxs.append(i + 1)

        qtine_slices = []
        starts = [0] + quarantine_idxs
        ends = quarantine_idxs + [len(self.I)]
        qtine_slices = list(zip(starts, ends))

        slice_widths = []
        # And then get max and width of each slice 
        for q_slice in qtine_slices:
            max_pair = max(enumerate(self.I[q_slice[0]:q_slice[1]]), 
                           key=lambda p: p[1])
            half_max = max_pair[1] / 2
            half_max_idxs = []
            for i in range(q_slice[0], q_slice[1] -1 ):
                if (self.I[i] -  half_max) * (self.I[i+1] - half_max) <= 0:
                    half_max_idxs.append(i)
            lo_t = self.t[half_max_idxs[0]]
            hi_t = self.t[half_max_idxs[-1]]
            slice_widths.append(hi_t - lo_t)
        return slice_widths



class AggregateTuple:
    def __init__(self, tup_list):
        """ Placeholder class to hold/access lists of TupleSIR objects """
        self.tup_list = tup_list 

    @classmethod
    def from_binary(cls, binary_obj):
        return cls(pickle.loads(binary_obj))

    def to_binary(self):
        return Binary(pickle.dumps([_.to_dict() for _ in self.tup_list]))

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
        


    def gather_agg(self, plot_series='I'):
        """ Returns a list of [(t, mean_I, std_I)] 
            for either I or R series, where the T's are from the one that has 
            the max (noninf) time 
        """
        def getter(tup):
            if plot_series == 'I':
                return self[idx].I
            else:
                return self[idx].R
        max_t_idx = max(enumerate(self), key=lambda tup: tup[1].t[-2])[0]
        ts = self[max_t_idx].t
        # First collect into lists 
        output = [] 
        for t in ts: 
            sublist = []
            for tup in self:
                els = tup.I if (plot_series == 'I') else tup.R 
                sublist.append(utils.linear_interp(tup.t, els, t)) 
            output.append(sublist) 
        # Hm... maybe better to return the full list and then decide from there
        return [(ts[i], _) for i, _ in enumerate(output)]
        #return [(ts[i], np.mean(_), np.std(_)) for i, _ in enumerate(output)]





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
    return G, TupleSIR.from_summary(summary, tmax=tmax), summary


def run_until_prop_IR(G, tau, gamma, rho, tmax, prop, total_nodes=None,
                      term_I=False, return_summary=True):
    """
    Runs SIR model until prop (in [0,1]) fraction of nodes are in I+R states 
    If total_nodes is not None, then the proportion is WRT total_nodes (and not len(G))
    
    RETURNS: G2, status_dict
    G2 is a copy of the graph G with I,R nodes removed 
    status_dict is a dict with keys {t, S, I, R} pointing to the right times 
    """

    if min([tmax, prop]) == 0:
        if return_summary:
            empty_summary = EoN.fast_SIR(G, tau, gamma, rho=0, tmax=0, return_full_data=True)

            return G, TupleSIR.empty(), empty_summary
        else:
            return G, TupleSIR.empty()

    total_nodes = total_nodes or len(G)
    re_init_count = rho * total_nodes / len(G)
    threshold = total_nodes * prop
    kwarg = {'term_IR': threshold}    
    if term_I:
        kwarg = {'term_I': threshold}
    # This has to be slower because we need to run the infection 
    # and then figure out which time to cut things off (ex-post facto)
    try:
        # IF USING MATT'S HACKED/FASTER VERSION OF EoN
        summary = EoN.fast_SIR(G, tau, gamma, rho=re_init_count, tmax=tmax, return_full_data=True,
                               **kwarg)
        breaktime = summary.t()[-1]
        summary_dict = utils.invert_dict(summary.get_statuses(time=breaktime))        
        G = G.copy() 
        G.remove_nodes_from(summary_dict.get('I', []))
        G.remove_nodes_from(summary_dict.get('R', []))
        if return_summary:
            return G, TupleSIR.from_summary(summary), summary
        else:
            return G, TupleSIR.from_summary(summary)
    except Exception as err:
        print("EXCEPTION", err)
        print("^Probably not using the local EoN")


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
    print("TRUNCI", trunc_I)
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


def quarantine_by_prop(G, tau, gamma, rho, prop_list, tmax, num_iter=1,
                       term_I=False, return_summary=False):
    # SINGLE QUARANTINE ONLY!!!
    # Runs a single quarantine which is initialized when the proportion I,R gets to prop
    # and then runs until tmax afterwards

    if not isinstance(prop_list, list):
        prop_list = [prop_list]

    outputs = [] 
    summaries = []
    original_G = G
    init_infect = round(rho * len(G))
    new_rho = init_infect / len(G)
    cured = False
    for iter_num in range(num_iter):
        count = 0
        G = original_G 
        tups = [] 
        remaining_time = tmax
        for prop in prop_list:
            if prop == 0:
                continue

            prop_out = run_until_prop_IR(G, tau, gamma, new_rho, remaining_time, 
                                       prop, total_nodes=len(original_G),
                                       term_I=term_I, return_summary=return_summary)
            if not return_summary:
                G, tup = prop_out 
            else:
                G, tup, summary = prop_out
                summaries.append(summary)
            new_rho = init_infect / len(G)
            count += 1
            remaining_time -= tup.tmax
            tups.append(tup) 
            if tup.I[-1] == 0: 
                cured = True
                break

        if remaining_time > 0.0 and not cured:
            run_until_time_out = run_until_time(G, tau, gamma, new_rho, remaining_time)
            tups.append(run_until_time_out[1])
            summaries.append(run_until_time_out[-1])
            count +=1
        outputs.append(TupleSIR.cat(tups))
        outputs[-1].count = count

    if len(outputs) == 1:
        if return_summary:
            return outputs[0], summaries
        return outputs[0]
    else:
        if return_summary:  
            return AggregateTuple(outputs), summaries
        return AggregateTuple(outputs)


# ======================================================================
# =           PLOTTING METHODS                                         =
# ======================================================================




SERIES_IDX = {'S': 1, 'I': 2, 'R': 3, 'Z': 2}
def plot_vanilla_run(G, tau, gamma, rho, tmax, series='IR', num_iter=5):
    axlist = []
    for serie in series: 
        fig, ax = plt.subplots(figsize=(8,8))
        for i in range(num_iter):
            runtup = EoN.fast_SIR(G, tau, gamma, rho=rho, tmax=tmax)
            if serie == 'Z':
                ax.plot(runtup[0], [sum(_)for _ in zip(runtup[SERIES_IDX['I']],
                                                       runtup[SERIES_IDX['R']])],
                        c='k', alpha=0.3)
            else:
                ax.plot(runtup[0], runtup[SERIES_IDX[serie]], c='k', alpha=0.3)
        ax.set_xlabel('Time')
        ax.set_ylabel('# %s' % serie)
        ax.set_title('%s vs time' % serie)
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
        print("Quarantine @ props p=%s,%s" % (p0, p1))
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



def heatmapify(grid, grid_idxs, axlabels='time', title=None):
    """ Builds a heatmap and formats it (except for the title) """
    assert axlabels in ['time', 'prop']
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.heatmap(grid, ax=ax)
    yticks = [_[0][0] for _ in grid_idxs]
    xticks = [_[1] for _ in grid_idxs[0]]
    if axlabels == 'time':
        ax.set_ylabel("Time of first quarantine")
        ax.set_xlabel("Time of second quarantine (after first)")
    else:
        ax.set_ylabel('Proportion of IR at Q1')
        ax.set_xlabel('Proportion of IR at Q2')
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
class DegreeXDict:
    def __init__(self, degree_x):
        self.degree_x = degree_x

    @classmethod
    def from_start_end(cls, start_G, end_G):
        """ Method to get data dicts: 
        ARGS:   
            start_G: graph before running SIR-quarantine 
            end_G: graph after running SIR-quarantine 
        RETURNS:         
            dict like {d: {'original': num,
                           'final': num,
                           N_i: num}, ...}
        Where d := degree at start-graph
               'original' -> number of nodes with degree d at start 
               'final' -> number of nodes that started with degree d 
                          that haven't been infected by end
                N_i -> int that points to how many nodes of degree d end up 
                       with degree d_i 
        """
        # Helper setups
        invert_start = utils.invert_dict(dict(start_G.degree)) # deg -> list of nodes
        end_degree = dict(end_G.degree) # node -> degree of that node 

        # Make output dict and loop over start-degrees
        output_dict = {}

        for start_deg, node_list in invert_start.items():
            saved_nodes = {node: end_degree[node] for node in node_list
                           if node in end_degree}
            key_dict = {'original': len(node_list),
                        'final': len(saved_nodes)}
            invert_end = {k: len(v) for k,v in utils.invert_dict(saved_nodes).items()}
            key_dict.update(invert_end)
            output_dict[start_deg] = key_dict

        return cls(output_dict)

    @classmethod
    def from_start_ends(cls, start_G, end_Gs):
        instances = [cls.from_start_end(start_G, _) for _ in end_Gs]
        return cls.merge(instances)


    @classmethod 
    def merge(cls, instances):
        return cls(utils.mergesum([_.degree_x for _ in instances]))


    def percent_survived(self):
        """ Returns a dict mapping start_degree-> %-susceptible after running
            For each dict, this is final / original
        """

        op = lambda d: d['final'] / d['original']
        return {k: op(v) for k,v in self.degree_x.items()}


    def avg_degrees(self):
        def op(v):
            # Computes average degree of nodes that remain 
            # (nodes that are removed (aka in IR) are not counted here
            num_nodes = v['final']
            if num_nodes == 0:
                return 0
            runsum = 0
            for key, subv in v.items():
                if key in ['final', 'original']:
                    continue 
                runsum += subv * key
            return runsum / num_nodes

        return {k: op(v) for k,v in self.degree_x.items()}

    def avg_degree_drop_percent(self):
        return {k: 1 - v /k for k,v in self.avg_degrees().items()}



