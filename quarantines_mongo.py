"""
General helpers for writing/reading experiment data from mongodb 
(because I like some persistence, it'll help aggregate data for me)


BASIC MONGODB COLLECTION SCHEMA:
DATABASE : quarantine20
	graphs: 
		process : str - like Barabasi_albert, etc etc 
		N 		: int - number of nodes
		params  : object - all args that go into graph construction
		seed : int - random seed used to generate this graph
	
	epidemics: 
		tau : float - infection rate
		gamma : float - recovery rate 
		rho:    float - params for initializing an infection
		name : str - shorthand helper?

	quarantine_by_props:
		graph_id : id of graph used 
		epidemic_id : id of epidemic used 
		quarantine_props : float[] - floats of when the quarantine was run 
		iter_num : int - which iteration this is 
		final_R: float - final number of recovered individuals 
		max_I : float - maximum number of infected individuals 
		full_data : object - {t, S, I, R} lists, keyed by tSIR

	quarantine_by_times: 
		graph_id : id of graph used 
		epidemic_id : id of epidemic used 
		quarantine_times : float[] - floats of when the quarantine was run 
		iter_num : int - which iteration this is 
		final_R: float - final number of recovered individuals 
		max_I : float - maximum number of infected individuals 
		full_data : object - {t, S, I, R} lists, keyed by 

	degree_obj: 
		graph_id : id of graph used 
		epidemic_id : id of epidemic used 


	vanilla_runs:
		graph_id: id of graph used 
		epidemic_id: id of epidemic used 
		num_runs: int - number of runs 
		final_Rs: int[] - list of ints with final recovered numbers
		max_Is: int[] - list of ints with maximum number of infected nodes 
		full_data: AggregateTuple - pickled aggregateTuple method


	plc_qtines:
		graph_id: id of graph used 
		epidemic_id: id of epidemic used 
		num_runs: int - number of trials 
		num_qtines: int - number of quarantines performed 
		qtine_strat: str - either 'prop' or 'time'
		qtine_data: float[] - list of when quarantines were applied 
		final_Rs: int[] - list of ints with final recovered numbers 
		max_Is: int[] - list of ints with max infected numbers

"""


from pymongo import MongoClient
from bson.objectid import ObjectId
from bson.binary import Binary
import quarantines as Q 
import networkx as nx
import random
import pickle
DATABASE_NAME = 'quarantine20'

# HELPERS: 

def insert_ba_graph(client, N, m):
	collection = client[DATABASE_NAME].graphs 
	seed = random.randint(1, 2 ** 20)
	doc = {'process': 'barabasi_albert', 
		   'N': N, 
		   'params': {'m': m}, 
		   'seed': seed}
	return collection.insert(doc)

def insert_plc_graph(client, N, m, p):
	collection = client[DATABASE_NAME].graphs 
	seed = random.randint(1, 2 ** 20) 
	doc = {'process': 'powerlaw_cluster', 
	       'N': N, 
	       'params': {'m': m, 'p': p},
	       'seed': seed}	       
	return collection.insert(doc)


def recreate_graph(client, graph_id):
	doc = client[DATABASE_NAME].graphs.find_one(ObjectId(graph_id))
	if doc['process'] == 'barabasi_albert':
		graph = nx.barabasi_albert_graph(doc['N'], doc['params']['m'], 
										 seed=doc['seed'])
	elif doc['process'] == 'powerlaw_cluster':
		graph = nx.powerlaw_cluster_graph(doc['N'], doc['params']['m'], 
										  doc['params']['p'], 
										  seed=doc['seed'])
	else:
		raise NotImplementedError

	graph.graph_id = doc['_id']
	return graph


def insert_epidemic_params(client, tau, gamma, rho, name=None):
	collection = client[DATABASE_NAME].epidemics
	doc = {'tau': tau, 'gamma': gamma, 'rho': rho}

	# First check to see if epidemic params exist: 
	checkdoc = collection.find_one({'tau': tau, 'gamma': gamma, 'rho': rho})
	if checkdoc is not None:
		return checkdoc['_id']

	if name is not None:
		doc['name'] = name
	return collection.insert(doc)

def collect_epidemic_params(client, epidemic_id):
	doc = client[DATABASE_NAME].epidemics.find_one(ObjectId(epidemic_id))
	return {k: doc[k] for k in ['tau', 'gamma', 'rho'] if k in doc}


# Populate quarantines:

def quarantine_by_prop_doc(graph, epidemic_id, epidemic_params, prop_list,
						   which_iter):
	# Get the python output objects
	tup = Q.quarantine_by_prop(graph, **epidemic_params, tmax=float('inf'),
 						       prop_list=prop_list, num_iter=1)

	# Process into mongo insert objects 
	output_doc = {'graph_id': graph.graph_id,
				  'epidemic_id': epidemic_id,
				  'quarantine_props': prop_list, 
				  'iter_num': which_iter, 
				  'final_R': tup.get_final_R().item(), 
				  'max_I': tup.get_max_I().item(), 
				  'full_data': Binary(pickle.dumps(tup.to_dict()))}

	return output_doc


def quarantine_by_time_doc(graph, epidemic_id, epidemic_params, time_list, 
							which_iter):
	tup = Q.quarantines_by_time(graph, **epidemic_params,
							    qtimes=time_list, tmax=float('inf'), num_iter=1)

	output_doc = {'graph_id': graph.graph_id,
				  'epidemic_id': epidemic_id,
				  'quarantine_times': time_list, 
				  'iter_num': which_iter, 
				  'final_R': tup.get_final_R().item(), 
				  'max_I': tup.get_max_I().item(), 
				  'full_data': Binary(pickle.dumps(tup.to_dict()))}
	return output_doc


def populate_quarantine_by_prop(client, graph_id, epidemic_id, 
								prop_range, iter_num):
	graph = recreate_graph(client, graph_id)
	epidemic_params = collect_epidemic_params(client, epidemic_id)
	collection = client[DATABASE_NAME]['quarantine_by_props']
	for prop_list in prop_range:
		for which_iter in range(iter_num):
			doc = quarantine_by_prop_doc(graph, epidemic_id, 
										 epidemic_params, prop_list,
									     which_iter)
			collection.insert_one(doc)

def populate_quarantine_by_time(client, graph_id, epidemic_id, time_range, 
								iter_num):
	graph = recreate_graph(client, graph_id)
	epidemic_params = collect_epidemic_params(client, epidemic_id)
	collection = client[DATABASE_NAME]['quarantine_by_times']
	for time_list in time_range:
		for which_iter in range(iter_num):
			doc = quarantine_by_time_doc(graph, epidemic_id, 
										 epidemic_params, time_list,
									     which_iter)
			collection.insert_one(doc)




def populate_percent_survived_by_time(client, graph_id, epidemic_id, 
									  time_ranges=None, prop_ranges=None, 
									  iter_num=100):
	""" Okay, now the real reason why I built this stupid mongo setup: 
		Function that runs |iter_num| simulations for each time in time ranges 
		where each simulation builds the DegreeXDict.percent_survived(...) 
		results... But I guess I want the following schema:

		Ultimate end goal is to plot the following (for various graphs/epidemics):
		- fix runtime, plot (start_degree, %alive-at-end)
		- memory issues here, so want documents like: 
			+ graph_id, epidemic_id, run_id, stop_time/stop_prop
			+ degree_start + degree_end 
			(maybe one doc per run) [ can alawys handle one doc in memory]

	ARGS:
		client: mongo client 
		graph_id : id of the graph stored in mongo 
		epidemic_id: id of the epidemic parameters stored in mongo 

	"""
	run_id = ObjectId()
	G = recreate_graph(client, graph_id)
	epidemic_params = collect_epidemic_params(client, epidemic_id)
	run_id = ObjectId()	
	assert (time_ranges is None) + (prop_ranges is None) == 1


	def doc_maker(G, stop_time=None, stop_prop=None, 
				  graph_id=graph_id, epidemic_id=epidemic_id, 
				  epidemic_params=epidemic_params, run_id=run_id,
				  which_iter=None):
		# makes a single document to insert into the collection
		base_dict = {'graph_id': graph_id, 
					 'epidemic_id': epidemic_id, 
					 'run_id': run_id, 
					 'iter': which_iter
					}

		if stop_time is not None:
			out_G = Q.run_until_time(G, **epidemic_params, tmax=stop_time)[0]
			base_dict['stop_time'] = stop_time[0]
		else:
			assert stop_prop is not None 
			out_G = Q.run_until_prop_IR(G, **epidemic_params, tmax=float('inf'), 
										prop=stop_prop)[0]
			base_dict['stop_prop'] = stop_prop
		degree_dict = Q.DegreeXDict.from_start_end(G, out_G)
		out_dict = {} 
		for k, v in degree_dict.degree_x.items():
			out_dict[str(k)] = {k2: v[k2] for k2 in ['original', 'final']}
		base_dict['survived'] = out_dict 
		print(base_dict)
		return base_dict

	collection = client[DATABASE_NAME].percent_survived
	for which_iter in range(iter_num):
		if time_ranges is not None:
			for t in time_ranges:
				doc = doc_maker(G, stop_time=t, which_iter=which_iter)
				collection.insert_one(doc)
		if prop_ranges is not None:
			for p in prop_ranges:
				doc = doc_maker(G, stop_prop=p, which_iter=which_iter)
				collection.insert_one(doc)



def populate_vanilla_runs(client, graph_id, epidemic_id, num_runs):
	graph = recreate_graph(client, graph_id)
	epidemic_params = collect_epidemic_params(client, epidemic_id)
	collection = client[DATABASE_NAME]['vanilla_runs']
	runs = [Q.run_until_time(graph, **epidemic_params, tmax=float('inf'))[1] 
			for _ in range(num_runs)]

	max_Is = [_.get_max_I().item() for _ in runs]
	final_Rs = [_.get_final_R().item() for _ in runs]
	agg_tup = Q.AggregateTuple(runs)

	doc = {'graph_id': ObjectId(graph_id), 
	   	   'epidemic_id': ObjectId(epidemic_id),
	   	   'num_runs': num_runs, 
	   	   'max_Is': max_Is,
	   	   'final_Rs': final_Rs, 
	   	   'full_data': agg_tup.to_binary()}
	return collection.insert(doc)


def populate_plc_qtines(client, graph_id, epidemic_id, num_runs, qtine_props):
	if not isinstance(qtine_props, list):
		qtine_props = [qtine_props]

	graph = recreate_graph(client, graph_id)
	epidemic_params = collect_epidemic_params(client, epidemic_id)
	collection = client[DATABASE_NAME]['plc_qtines']
	agg_tup = Q.quarantine_by_prop(graph, **epidemic_params, prop_list=qtine_props, 
								   tmax=float('inf'), num_iter=num_runs)
	max_Is = [_.get_max_I().item() for _ in agg_tup.tup_list]
	final_Rs = [_.get_final_R().item() for _ in agg_tup.tup_list]
	doc = {'graph_id': ObjectId(graph_id), 
		   'epidemic_id': ObjectId(epidemic_id), 
		   'num_runs': num_runs, 
		   'num_qtines': len(qtine_props), 
		   'qtine_strat': 'prop',
		   'qtine_data': qtine_props, 
		   'final_Rs': final_Rs, 
		   'max_Is': max_Is}
	return collection.insert(doc)



if __name__ == '__main__':
	main()

