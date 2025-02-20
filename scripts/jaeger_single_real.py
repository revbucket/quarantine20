import sys 
sys.path.append('..')
import networkx as nx
import EoNlocal as EoN
import quarantines as Q 
import utils
import random
import matplotlib.pyplot as plt
import gc
import graph_generators as gg
import quarantines_mongo as qm 
from pymongo import MongoClient


import argparse 

RESOLUTION = 100

def main(args):
	name = args.name
	# Parse graph params
	graph_type = args.graph_type 
	graph_params = args.graph_params
	if len(graph_params) == 1:
		graph_params = graph_params[0]
	#Parse infection params
	tau = float(args.tau )
	gamma = float(args.gamma )



	# Parse run params 
	num_iter = args.num_iter

	# Run thing we want to run (and print as we go along!) 
	db = MongoClient('localhost', 27017)[qm.DB2]
	if graph_type == 'gemsec_fb':
		graph = gg.load_gemsec_fb(graph_params)

	elif graph_type == 'gemsec_deezer':
		graph = gg.load_gemsec_deezer(graph_params)

	elif graph_type == 'arxiv':
		graph = gg.load_arxiv_collab(graph_params)

	elif graph_type == 'hiv':
		graph = gg.load_hiv(int(graph_params))
		return graph
	elif graph_type == 'hs':
		minutes = float(graph_params[0])
		scale = int(graph_params[1])
		graph = gg.load_highschool(minutes, scale=scale)
		return graph

	else:
		raise NotImplementedError("Nope")
	rho = float(args.init_infected / len(graph))
	epi_id = qm.insert_epidemic_params(db, tau, gamma, rho, name=name)


	prop_range = [_ / RESOLUTION for _ in range(0, RESOLUTION)]

	qm.populate_quarantine_by_prop(db, None, epi_id, prop_range, num_iter, 
									   name=name, save_full_data=False, graph=graph)




if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Process some integers.')

	parser.add_argument('-g', '--graph_type', dest='graph_type', type=str, 
						choices=['gemsec_fb', 'gemsec_deezer', 'arxiv', 'hiv', 'hs'], 
		                required=True)
	parser.add_argument('-gparam', '--graph_params', dest='graph_params', type=str, 
					     							 nargs='+', required=True)
	parser.add_argument('-tau', '--tau', dest='tau', type=float, required=False, default=0.1)
	parser.add_argument('-gamma', '--gamma', dest='gamma', type=float, required=False, default=0.1)
	parser.add_argument('-init', '--init_infected', dest='init_infected', type=int, required=False, default=10)
	parser.add_argument('-iter', '--num_iter', dest='num_iter', type=int, required=False, default=20)
	parser.add_argument('-name', '--name', dest='name', type=str, required=True)


	args = parser.parse_args()	
	output = main(args)

