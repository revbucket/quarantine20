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

def ba10(name, m, N, tau, gamma, num_init, num_iter):
	print("RUNNING BA", locals())
	db = MongoClient('localhost', 27017)[qm.DB3] 
	rho = float(num_init / N)
	graph_id = qm.insert_ba_graph(db, N, m, name=name)
	epi_id = qm.insert_epidemic_params(db, tau, gamma, rho, name=name)


	prop_list = []
	for i in range(RESOLUTION):
		for j in range(RESOLUTION-i):
			prop_list.append([i / RESOLUTION, j / RESOLUTION])
	print("NEUTERED")
	return 
	qm.populate_quarantine_by_prop(db, graph_id, epi_id, prop_list, num_iter, 
								   name=name, save_full_data=False)


def fb_artist(name, tau, gamma, num_init, num_iter):
	print("RUNNING GEMSEC ARTIST", locals())
	db = MongoClient('localhost', 27017)[qm.DB3] 

	graph = gg.load_gemsec_fb('artist')	
	rho = float(num_init / len(graph))
	epi_id = qm.insert_epidemic_params(db, tau, gamma, rho, name=name)
	

	prop_list = []
	for i in range(RESOLUTION):
		for j in range(RESOLUTION-i):
			prop_list.append([i / RESOLUTION, j / RESOLUTION])
	print("NEUTERED")
	return 
	qm.populate_quarantine_by_prop(db, None, epi_id, prop_list, num_iter, 
								   name=name, save_full_data=False, graph=graph)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-g', dest='graph_type', type=str, choices=['fb.artist', 'ba10'])

	args = parser.parse_args() 

	if args.graph_type == 'fb.artist':
		fb_artist('fb_artist_2Q', 0.5, 1.0, 10, 25)
	if args.graph_type == 'ba10':
		ba10('ba10_2Q', 10, 1e4, 0.5, 1.0, 10, 25)