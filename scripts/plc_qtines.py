import sys 
sys.path.append('..')

import pymongo 
import networkx as nx 
import numpy as np 
import matplotlib.pyplot as plt 

import quarantines as Q 
import quarantines_mongo as qm 

DATABASE_NAME = qm.DATABASE_NAME 



def prod(iter_list):
	count = 1 
	for el in iter_list:
		count *= len(el)
	return count 

if __name__ == '__main__':
	client = pymongo.MongoClient('localhost', 27017)
	N = 10 ** 4
	ms = [3, 5] 
	p = 0.4
	tau = 0.1 
	r0s = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0] 
	rho = 25 / N
	qprop = [[x / 20, y / 20] for x in range(16) for y in range(16)]
	num_runs = 5 

	total_count = prod([ms, r0s, qprop])
	count = 0 
	for m in ms:
		graph_id = qm.insert_plc_graph(client, N, m, p)
		for r0 in r0s:
			epidemic_id = qm.insert_epidemic_params(client, tau, tau / r0, rho)
			for qlist in qprop:
				count += 1
				print("%s/ %s" % (count, total_count))
				qm.populate_plc_qtines(client, graph_id, epidemic_id, num_runs, 
									   qlist)
	main()