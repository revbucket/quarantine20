""" Python script that should be run ON Jaeger -- does the following: 
	- Builds a graph and epidemic parameters according to commandline args 
	- Runs populate_quarantine_by_prop for proprange [0.01 * _ for _ in range 0, 100]
	  and iter_num specified by command line 
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


import argparse 

RESOLUTION = 10

def main(args):
	name = args.name
	# Parse graph params
	graph_type = args.graph_type 
	N = int(float(args.N))
	graph_params = args.graph_params 
	
	# Parse infection params  		
	tau = float(args.tau )
	gamma = float(args.gamma )
	rho = float(args.init_infected / N)


	# Parse run params 
	num_iter = args.num_iter

	# Run thing we want to run (and print as we go along!) 
	db = MongoClient('localhost', 27017)[qm.DB2]
	if graph_type == 'ba':
		m = int(graph_params[0])
		graph_id = qm.insert_ba_graph(db, N, m, name=name)
	elif graph_type == 'plc':
		m = int(graph_params[0])
		p = float(graph_params[1])
		graph_id = qm.insert_plc_graph(db, N, m, p, name=name)
	elif graph_type == 'rw':
		qe = float(graph_params[0])
		qv = float(graph_params[1])
		graph_id = qm.insert_rw_graph(db, N, qe, qv, name=name)
	elif graph_type == 'nn':
		u = float(graph_params[0])
		k = int(graph_params[1])
		graph_id = qm.insert_nn_graph(db, N, u, k, name=name)
	elif graph_type == 'ws':
		k = int(graph_params[0])
		p = float(graph_params[1])
		graph_id = qm.insert_ws_graph(db, N, k, p, name=name)
	else:
		raise NotImplementedError("Only BA/PLC graphs supported")

	epi_id = qm.insert_epidemic_params(db, tau, gamma, rho, name=name)


	prop_range = [_ / RESOLUTION for _ in range(0, RESOLUTION)]

	qm.populate_quarantine_by_prop(db, graph_id, epi_id, prop_range, num_iter, 
									   name=name, save_full_data=False)




if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Process some integers.')

	parser.add_argument('-g', '--graph_type', dest='graph_type', type=str, 
						choices=['ba', 'plc', 'rw', 'nn', 'ws'], 
		                required=True)
	parser.add_argument('-N', '--N', dest='N', required=True)
	parser.add_argument('-gparam', '--graph_params', dest='graph_params', nargs='+', required=True)
	parser.add_argument('-tau', '--tau', dest='tau', type=float, required=True)
	parser.add_argument('-gamma', '--gamma', dest='gamma', type=float, required=True)
	parser.add_argument('-init', '--init_infected', dest='init_infected', type=int, required=True)
	parser.add_argument('-iter', '--num_iter', dest='num_iter', type=int, required=True)
	parser.add_argument('-name', '--name', dest='name', type=str, required=True)


	args = parser.parse_args()	
	main(args)

