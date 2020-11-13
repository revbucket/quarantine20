import sys
import quarantines as Q
import quarantines_mongo as qm
import graph_generators as gg 
import networkx as nx 
import networkit as nk
import pymongo
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from prettytable import PrettyTable
from tabulate import tabulate
from scipy import stats
from pprint import pprint
import os 
import random
sns.set()
DIRNAME = os.path.dirname(__file__)
DATA_NAME = 'quarantine_by_props_1027.pkl'


def gather_data(data_name=DATA_NAME):
    data_file = open(os.path.join(DIRNAME, data_name), 'rb')
    data = pickle.load(data_file)
    data_by_name = {}
    for datum in data:
        name = datum.get('name')
        if name not in data_by_name:
            data_by_name[name] = [] 
        data_by_name[name].append(datum)
    return data_by_name
data = gather_data()


# Part 0: Collect data by name and print all names...

def group_by_qprop(series):
    gather = {} 
    for datum in series:
        qprop = datum['quarantine_props']
        if qprop not in gather:
            gather[qprop] = [] 
        gather[qprop].append(datum)
        
    # and now modify to only collect final_R, max_I
    output = {}
    for qprop, datalist in gather.items():
        final_rs = [_['final_R'] for _ in datalist]
        max_is = [_['max_I'] for _ in datalist]
        output[qprop] = {'final_R': final_rs, 
                         'max_I': max_is}
    return output

def gather_by_name(data, name, i_or_r):
    # Gather triples of [(qprop, mean, std), ...]
    assert i_or_r in ['I', 'R']
    series = group_by_qprop(data[name])
    output = []
    def getter(doc):
        if i_or_r == 'I':
            return doc['max_I']
        else:
            return doc['final_R']
        
        
    for k, v in series.items():
        mean = np.mean(getter(v))
        std = np.std(getter(v))
        output.append((k, mean, std))
        
    return sorted(output, key=lambda trip: trip[0])
    

def size_lookup_by_name(name):
    if name in SYNTHETICS:
        return 10 ** 4 
    elif name.startswith('fb.'):
        return len(gg.load_gemsec_fb(name.split('.')[1]))
    elif name.startswith('deezer.'):
        return len(gg.load_gemsec_deezer(name.split('.')[1]))
    elif name.startswith('arxiv.'):
        return len(gg.load_arxiv_collab(name.split('.')[1]))
    elif name.startswith('hiv') or name.startswith('hs'):
        G = recreate_by_name(name)
        return len(G)
    
    
def recreate_by_name(name):
    # Don't worry too much about actual random seed, just gather parameters
    if name.startswith('ba'):
        if name.startswith('ba10'):
            return gg.ba_graph(10 **4, 10)
        return gg.ba_graph(10 ** 4, int(name[2:]))
    elif name.startswith('plc'):
        m, p = name[3:].split('.')
        return gg.plc_graph(10 **4, int(m), float(p) / 100)
    elif name.startswith('rw'):
        assert name == 'rw.91.94'
        return gg.random_walk_graph(10 ** 4, 0.91, 0.94)
    elif name.startswith('nn'):
        assert name == 'nn.886'
        return gg.nearestNeighbor_mod(10 **4, 0.88,6)
    elif name.startswith('ws'):
        assert name == 'ws10.05'
        return gg.watts_strogatz_graph(10 ** 4, 10, 0.05)
    elif name.startswith('fb.'):
        return gg.load_gemsec_fb(name.split('.')[1])
    elif name.startswith('deezer.'):
        return gg.load_gemsec_deezer(name.split('.')[1])
    elif name.startswith('arxiv.'):
        return gg.load_arxiv_collab(name.split('.')[1])
    elif name.startswith('hs'):
        params = name[2:].split('_')
        return gg.load_highschool(float(params[0]), int(params[1]))
    elif name.startswith('hiv'):
        param = int(name[3:].split('R')[0])
        return gg.load_hiv(param)
    elif name.startswith('ba10'):
        return gg.ba_graph(10 **4, 10)
    elif name.startswith('fbartist'):
        return gg.load_gemsec_fb('artist')
    elif name.startswith('gemsec_fb'):
        if 'ARTIST' in name:
            return gg.load_gemsec_fb('artist')

        return gg.load_gemsec_fb(name[10:-6])
    elif name.startswith('gemsec_deezer'):
        return gg.load_gemsec_deezer(name.split('_')[2])
    elif name.startswith('arxiv_'):
        return gg.load_arxiv_collab(name.split('_')[1])

    else:
        raise Exception("wut")
    return
    
    
def parse_data_name(name):
    prefixes = ['ba', 'plc', 'nn', 'rw', 'ws' ]
    match = None
    for prefix in prefixes:
        if name.startswith(prefix):
            match = prefix 
    if match is None or len(name.split('_')) == 1:
        return recreate_by_name(name)
    params = name.split('_')
    R = params[-1]
    if name.startswith('ba'):
        if name.startswith('ba10'):
            return gg.ba_graph(10 ** 4, 10)
        return gg.ba_graph(10 ** 4, int(params[1]))
    elif name.startswith('plc'):
        return gg.plc_graph(10 **4, int(params[1]), float(params[2]))
    elif name.startswith('rw'):
        return gg.random_walk_graph(10 ** 4, float(params[1]), float(params[2]))
    elif name.startswith('nn'):
        return gg.nearestNeighbor_mod(10 **4, float(params[1]),int(params[2]))
    elif name.startswith('ws'):
        return gg.watts_strogatz_graph(10 ** 4, int(params[1]), float(params[2]))

    
def get_r_by_name(name):
    prefixes = ['ba', 'plc', 'nn', 'rw', 'ws' ]
    match = None
    for prefix in prefixes:
        if name.startswith(prefix):
            match = prefix 
    if match is None or len(name.split('_')) == 1:
        return 1
    params = name.split('_')
    return float(params[-1])



def get_minR_graph(data, name):
    # Gathers the graph (after minR optimal quarantine) 
    # 1) Get minR prop:
    minprop = min(gather_by_name(data, name, 'R'), key=lambda trip: trip[1])[0]
    
    
    # 2) Recreate graph and rerun 
    R = get_r_by_name(name)
    G = parse_data_name(name)
    outG = Q.run_until_prop_IR(G, R, 1, 10 / len(G), float('inf'), minprop)[0]
    return outG


#fig, ax = plt.subplots(figsize=(10,10)) #<---- general axis maker
def plot_vs(data, names, irs, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,10))
    if not isinstance(names, list):
        names = [names]
    def plot_single_v(name, ir, ax=ax):
        try: 
            size = size_lookup[name]
        except:   
            size = 10 ** 4 # size_lookup[name]
        i_trips = gather_by_name(data, name, 'I')
        r_trips = gather_by_name(data, name, 'R')
        qprops = [_[0] for _ in i_trips]
        Is = [_[1] / size for _ in i_trips]
        qprops_r = [_[0] for _ in r_trips]
        Rs = [_[1] / size for _ in r_trips]
        if 'I' in ir:
            ax.plot(qprops, Is, label=name)
        if 'R' in ir:
            ax.plot(qprops, Rs, label=name)
    
    
    for name in names:
        plot_single_v(name, irs, ax=ax)
    ax.legend()
    

def degree_hist(G):
    pairs = {}
    for _, d in G.degree():
        pairs[d] = pairs.get(d, 0) + 1
    return pairs

def return_cbin(G):
    # Returns best-fit powerlaw exponent by using cumulative binning 
    items = sorted(degree_hist(G).items(), key=lambda p:-p[0])
    cdf = []
    runsum = 0
    for deg, num in items:
        runsum += num 
        cdf.append((deg, runsum))
    xform = [(np.log(_[0]), np.log(_[1] / len(G))) for _ in cdf]
    return xform 



def degree_hist(G):
    pairs = {}
    for _, d in G.degree():
        pairs[d] = pairs.get(d, 0) + 1
    return pairs


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
    return 1 -slope, r_value

        
def collect_graph_data_by_name(name, data=data):
    # Returns graph [name, nodes, average degree, cluster coefficient, powerlaw exponent, best_iTup, bestrTup] 
    G = recreate_by_name(name)
    nodes = len(G)
    avg_deg = sum(dict(G.degree()).values()) / len(G)
    cc = nx.average_clustering(G)
    powerlaw = get_powerlaw_exponent_cbin(G)
    no_q = gather_by_name(data, name, 'R')[0][1] /len(G)
    best_i = min(gather_by_name(data, name, 'I'), key=lambda tri: tri[1])
    best_i_prop = best_i[0]
    best_i_val = best_i[1] / nodes
    best_r = min(gather_by_name(data, name, 'R'), key=lambda tri: tri[1])
    best_r_prop = best_r[0]
    best_r_val = best_r[1] / nodes
    return [name, nodes, avg_deg, cc, powerlaw, no_q, best_r_prop, best_r_val]


def tabulate_by_name(names):
    headers=['name', 'nodes', 'avg_deg', 'cc', 'powerlaw', 'no_Q', 'best_r_prop', 'best_r_val']
    graph_data = [collect_graph_data_by_name(_) for _ in names]
    print(tabulate(graph_data, headers=headers, floatfmt='.2f'))
    return graph_data