""" Utilities for quarantine repo:
    - Generic helpers for iterables
    - Helpers for plotting things
"""
import networkx as nx 
from scipy import stats
import random
import numpy as np
import matplotlib.colors as mcolors

def selector(els, idxs):
    for idx in idxs:
        yield els[idx]

def tuple_filter(tup_iter, idx):
    return [tup[idx] for tup in tup_iter]

def invert_dict(d):
    new_dict = {}
    for k,v in d.items():
        if v in new_dict:
            new_dict[v].append(k)
        else:
            new_dict[v] = [k]
    return new_dict

def argmax(iterable):
    # Returns max_idx, max_value
    max_val, max_idx = -float('inf'), None
    for i, el in enumerate(iterable):
        if el > max_val:
            max_val = el
            max_idx = i
    return max_idx, max_val

def mean(iterable, lambda_=None):
    count, runsum = 0, 0.0 
    if lambda_ is None:
        lambda_ = lambda x: x 
    for el in iterable:
        count += 1
        runsum += lambda_(el)
    return runsum / float(count)
    

def mergesum(dicts):
    """ Given a list/iterable of (nested) dicts, will merge them together 
        where merge at the base level means summing values for shared keys 
    """
    def looper(iter_input, output=None):
        if all(isinstance(_, (float, int)) for _ in iter_input):
            return sum(iter_input)
        # Collect by shared keys:
        shared_keys = {}
        for el in iter_input:
            for k,v in el.items():
                if k not in shared_keys:
                    shared_keys[k] = []
                shared_keys[k].append(v)

        return {k: looper(v) for k, v in shared_keys.items()}
    return looper(dicts)

def binsearch(els, val):
    # Finds the largest index of els such that els[idx] < val
    def subroutine(start_idx, end_idx, els=els, val=val):
        if end_idx - start_idx <= 1:
            return start_idx
        intermed = int((start_idx + end_idx) / 2)
        
        if els[intermed] <= val:
            return subroutine(intermed, end_idx)
        else:
            return subroutine(start_idx, intermed)
    idx = subroutine(0, len(els) -1)
    next_idx = min([len(els) -1, idx + 1])
    if els[next_idx] <= val:
        return next_idx
    return idx
        
 

def linear_interp(xlist, ylist, xval):
    # Finds the linear interpolation of a AggregateTuple
    assert len(xlist) == len(ylist)
    if xval > xlist[-1]:
        return ylist[-1]
    if xval < xlist[0]:
        return xlist[0]
    
    idx = binsearch(xlist, xval)
    if xlist[idx] == xval:
        return ylist[idx]
    
    prop = (xval - xlist[idx]) / (xlist[idx + 1] - xlist[idx])
    return ylist[idx] + prop * (ylist[idx + 1] - ylist[idx])
    


### Plotting helpers


def c(i):
    return list(mcolors.TABLEAU_COLORS)[i]
    # return 'bgrcmyk'[i]

def select_mean(triplist):
    # input is a list of (xval, ymean, ystd)
    # output is just (xval, ymean)
    return ([_[0] for _ in triplist], [_[1] for _ in triplist])


def plotfill_trips(triplist):
    # input is a list of (xval, ymean, ystd)
    # output is (xval, ymean-yst, ymean+ystd)
    return ([_[0] for _ in triplist], 
            [_[1] - _[2] for _ in triplist], 
            [_[1] + _[2] for _ in triplist])

#### GRAPH UTILS 

def degree_hist(G):
    pairs = {}
    for _, d in G.degree():
        pairs[d] = pairs.get(d, 0) + 1
    return pairs


# Want to do scatter plot over 'shortest path length' and 'reduction in quarantine' 
def avg_deg(G):
    return 2 * G.number_of_edges() / len(G)


def get_powerlaw_exponent_cbin(G):
    # Returns best-fit powerlaw exponent by using cumulative binning 
    items = sorted(degree_hist(G).items(), key=lambda p:-p[0])
    cdf = []
    runsum = 0
    for deg, num in items:
        runsum += num 
        cdf.append((deg, runsum))
    xform = [(np.log(_[0]), np.log(_[1] / len(G))) for _ in cdf]
    #plt.scatter(*zip(*xform))
    slope, intercept, r_value, p_value, std_err = stats.linregress(*zip(*xform))
    return 1 -slope


def get_largest_cc(G):
    ccs = nx.connected_components(G)
    maxlen = 0
    max_ccset = None 
    for conn in ccs:
        if len(conn) > maxlen:
            maxlen = len(conn)
            max_ccset = conn
    return max_ccset


def avg_shortest_path(G, samples=10 * 1000):
    # Takes weighted average of connected component shortest paths? 
    # Randomly samples pairs of nodes in the largest connected component and computes average path length 
    ccset = get_largest_cc(G)
    sampleset = random.choices(list(ccset), k=samples * 2)
    len_sum = 0 
    count = 0
    for i in range(0, len(sampleset), 2):
        len_sum += nx.shortest_path_length(G, source=sampleset[i], target=sampleset[i+1])
        count += 1
    return len_sum / count