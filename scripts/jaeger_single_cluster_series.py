"""" Runs the double quarantine methods for 
	BA10 graph and FB.Artist graph 
"""

import sys 
sys.path.append('..')
import networkx as nx
import EoNlocal as EoN
import quarantines as Q 
import utils
import random
import matplotlib.pyplot as plt
import gc
import quarantines_mongo as qm 
from pymongo import MongoClient
import graph_generators as gg

import argparse 

RESOLUTION = 100

def make_plc_cluster(N, m, cluster):
    low_cluster = gg.plc_graph(N, m, 0)
    high_cluster = gg.plc_graph(N, m, 1.0)
    lo, lo_val = 0, nx.average_clustering(low_cluster)
    hi, hi_val = 1, nx.average_clustering(high_cluster)
    assert lo_val <= cluster <= hi_val
    cur, cur_val = lo, lo_val
    while abs(cur_val - cluster) > 0.001:
        cur = (lo + hi) / 2 
        G = gg.plc_graph(N, m, cur)
        cur_val = nx.average_clustering(G)
        if cluster > cur_val:
            lo, lo_val = cur, cur_val 
        else:
            hi, hi_val = cur, cur_val
    return cur


def cluster_series():
	SERIES = [0.05, 0.1, 0.15, 0.2, 0.25]
	N = 10 * 1000 
	M = 10
	tau = 0.5 
	gamma = 1.0 
	NUM_INIT = 10 
	NUM_ITER = 25 

	db = MongoClient('localhost', 27017)[qm.DB4]
	rho = float(NUM_INIT / N)


	for serie in SERIES:
		NAME = 'cluster_series_' + str(serie)
		p = make_plc_cluster(N, M, serie)
		graph_id = qm.insert_plc_graph(db, N, M, p)
		epi_id = qm.insert_epidemic_params(db, tau, gamma, rho, name=NAME)
		prop_list = [] 
		for i in range(RESOLUTION):
			prop_list.append(i / RESOLUTION)
		qm.populate_quarantine_by_prop(db, graph_id, epi_id, prop_list, NUM_ITER,
									   name=NAME, save_full_data=False, 
									   prefix=str(serie))


def size_series():
	SERIES = [5000, 10000, 50 * 1000, 100 * 1000]
	M = 10
	tau = 0.5 
	gamma = 1.0 
	NUM_INIT = 10 
	NUM_ITER = 25 

	db = MongoClient('localhost', 27017)[qm.DB4]



	for serie in SERIES:
		N = serie
		rho = float(NUM_INIT / N)		
		NAME = 'size_series_' + str(serie)

		graph_id = qm.insert_ba_graph(db, N, M)
		epi_id = qm.insert_epidemic_params(db, tau, gamma, rho, name=NAME)
		prop_list = [] 
		for i in range(RESOLUTION):
			prop_list.append(i / RESOLUTION)
		qm.populate_quarantine_by_prop(db, graph_id, epi_id, prop_list, NUM_ITER,
									   name=NAME, save_full_data=False, 
									   prefix=str(serie))




if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', dest='graph_type', type=str, choices=['cluster', 'size'])

	args = parser.parse_args() 

	if args.graph_type == 'cluster':
		cluster_series()
	if args.graph_type == 'size':
		size_series()