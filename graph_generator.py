""" File to handle some shorthand graph generation techniques 

TODOS: Find good default parameters here 
"""

import networkx as nx 



def scale_free_community(N, d, m, Q):
	""" Builds a scale-free graph with community structure as 
		described by Sah et al. 
	"""
	pass # WAITING ON EMAIL BACK FROM PRATHA 


def log_normal_config(N, mu, sigma):
	""" Builds a configuration model with the given graph size and 
		log-normal parameters
	"""
	samples = np.random.lognormal(mean=mu, sigma=sigma, size=N)
	# Maybe there's a better way to do this:
	return nx.expected_degree_graph(samples, selfloops=False)

	rounded_w = [int(round(_)) for _ in w]
	option_2 = nx.configuration_model(rounded_w)