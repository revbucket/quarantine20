""" File that builds the graphs with various parameters """


import networkx as nx
import networkit as nit
import random 
import math 
import os 
import bisect
DIRNAME = os.path.dirname(__file__)

# ================================================================================
# =           Networkx graph creation techniques                                 =
# ================================================================================

def ba_graph(N, m, seed=None):
    """ Good parameters here:
    GOOD PARAMETERS: [5, 10]
    """
    return nx.barabasi_albert_graph(N, m, seed=seed)


def plc_graph(N, m, p, seed=None):
    """ Good params ??? (probably like the random walk graph params)"""
    return nx.powerlaw_cluster_graph(N, m, p, seed=seed)


def random_walk_graph(N, qe, qv, seed=None):
    """ Good params:
        [0.90, 0.23], [0.91, 0.94], [0.93, 0.99], [0.96, 0.90], [0.93, 0.98]
    """
    G = nx.Graph()
    G.add_edge(0, 1)
    random.seed(seed)

    def random_neighbor(G, v):
        return random.choice(list(G.neighbors(v)))

    def random_walk_add(G, i):
        vlist = [random.randint(0, i - 1)]
        G.add_edge(i, vlist[-1])
        while random.random() < qe:
            vlist.append(random_neighbor(G, vlist[-1]))
        edges = [] 
        for v in set(vlist):
            if random.random() < qv:
                edges.append((i, v))
        G.add_edges_from(edges)


    for i in range(2, N):
        random_walk_add(G, i)
    return G


def watts_strogatz_graph(N, k, p, seed=None):
    return nx.watts_strogatz_graph(N, k, p, seed=seed)





def forestFire_mod_burnProcedure(G, node1, node2, myDict, p):
    """
    Helper function for forestFire_mod
    Handles burning procedure recursively for the modified Forest Fire model
    Recurive depth is handled by global variable 'limit'
    """
    # Ensure recursion does not go too deep
    global limit
    limit += 1

    if limit >= 5000:
        return

    # Generate a non-zero random floating point between 0 and 1
    y = 0
    while y == 0:
        y = random.random()

    # How many neighbors to "burn"
    x = (int)(math.ceil((math.log10(y) / math.log10(p)) - 1))
    burn = [] # Keep track of which neighbors to burn
    nbors = list(G.neighbors(node1))

    # No neighbors to burn
    if len(nbors) == 0:
        return G

    # If there are fewer neighbors than needed to burn, burn all
    elif len(nbors) <= x:
        for i in range(0, len(nbors)):
            if nbors[i] not in myDict:
                burn.append(nbors[i])
                myDict[nbors[i]] = nbors[i]

    # Choose the 'x' amount of neighbors to burn
    else:
        for i in range(0, x):
            a = random.randrange(0, len(nbors))
            b = 0
            for j in range(0, i):
                while nbors[a] == burn[j] or (nbors[a] in myDict):
                    a = random.randrange(0, len(nbors))
                    if nbors[a] in myDict:
                        b += 1
                    if (len(nbors) - b) < x:
                        break

                if (len(nbors) - b) < x:
                    break

            if (len(nbors) - b) < x:
                break

            burn.append(nbors[a])
            myDict[nbors[a]] = nbors[a]

    # Burn
    for i in range(0, len(burn)):
        if burn[i] != node2:
            G.add_edge(node2, burn[i])

    # Repeat recursively
    for i in range(0, len(burn)):
        forestFire_mod_burnProcedure(G, burn[i], node2, myDict, p)

def forestFire_mod(n, p, seed=None):
    """
    Generates a graph based on a modified Forest Fire model.

    This is a modified version of the Forest Fire model
    that creates undirected edges in the edge creation process.[1]

    Input:
        n = number of nodes in the graph (integer)
        p = "burn" rate (floating point)
    
    Output:
        nx.Graph()
    """
    # Prepare graph
    G = nx.Graph()
    random.seed(seed)
    nodeCounter = 0 # Tracks next available node ID

    # Keep adding nodes until we have 'n' nodes
    while nodeCounter < n:
        # Recursion limit
        global limit
        limit = 0
        target = nodeCounter
        G.add_node(target)
        nodeCounter = len(G) # Update next available nodeID

        if nodeCounter == 1:
            continue

        # Select a random node from graph that is not the current 'target' node
        randomNode = random.randrange(0, nodeCounter)

        myDict = dict()
        G.add_edge(randomNode, target)
        myDict[randomNode] = randomNode

        # Start burning
        forestFire_mod_burnProcedure(G, randomNode, target, myDict, p)

    return G


def nearestNeighbor_mod(n, u, k, seed=None):
    """
    Generates a graph based on a modified Nearest Neighbor  model.

    This is a modified version of the Nearest Neighbor model
    that creates an undirected graph with power-law exponent
    between 1.5 and 1.75, matching that of online social networks.
    This is done so that each time a new node is added, k random
    pairs of nodes in the connected component of the graph are connected. [1]

    Input:
        n = number of nodes in the graph (integer)
        u = probability that determines if a new node is added or if a pair of 2 hop neighbors is connected (floating point)
        k = each time a new node is added, k pairs of random nodes in the connected component are connected (integer)
    Output:
        nx.Graph()

    GOOD PARAMS: [0.53, 1], [0.90, 2], [0.88, 6], [0.88, 10], [0.90, 5]
    """
    random.seed(seed)
    G = nx.Graph()
    nodeCounter = 0 # Keeps track of ID of next available node

    degreeArray = [0 for i in range(0,n)]
    d = []
    N = [0 for i in range(0,2)]

    while nodeCounter < n: # Until we reach 'n' nodes...
        if random.random() < u:
            if len(d) == 0:
                continue

            x = random.choice(d) # Pick a node from list
            N = random.sample(list(G.neighbors(x)), 2) # Pick 2 unique nodes in the list of neighbors

            if not G.has_edge(N[0], N[1]): # If no edge exists between the 2, connect
                G.add_edge(N[0], N[1])
                degreeArray[N[0]] += 1
                degreeArray[N[1]] += 1
                if degreeArray[N[0]] == 2:
                    bisect.insort(d, N[0])

                if degreeArray[N[1]] == 2:
                    bisect.insort(d, N[1])

        else:
            nodesSoFar = nodeCounter
            G.add_node(nodeCounter)
            nodeCounter += 1

            if nodeCounter == 1: # No use in continuing if there is only one node in the graph
                continue

            a = random.randrange(0, nodesSoFar) # Pick a node in the graph
            G.add_edge(a, nodesSoFar)
            degreeArray[a] += 1
            degreeArray[nodesSoFar] += 1

            if degreeArray[a] == 2:
                bisect.insort(d, a)
            if degreeArray[nodesSoFar] == 2:
                bisect.insort(d, nodesSoFar)

            for i in range(0, k): # Connect k random pairs in the graph
                N[0] = random.randint(0, nodeCounter - 1)
                N[1] = N[0]

                while N[1] == N[0]: # Ensure the two nodes are different (no self-loops)
                    N[1] = random.randint(0, nodeCounter - 1)

                if not G.has_edge(N[0], N[1]):
                    G.add_edge(N[0], N[1])

                    degreeArray[N[0]] += 1
                    degreeArray[N[1]] += 1

                    if degreeArray[N[0]] == 2:
                        bisect.insort(d, N[0])
                    if degreeArray[N[1]] == 2:
                        bisect.insort(d, N[1])

    return G



# =================================================================
# =           REAL GRAPH LOADERS                                  =
# =================================================================

GEMSEC_FB_ARGS = ['artist', 'athletes', 'company', 'government', 
                  'new_sites', 'politician', 'public_figure', 'tvshow']
def load_gemsec_fb(s):
    assert s in GEMSEC_FB_ARGS
    return nx.read_adjlist(os.path.join(DIRNAME, 'snap/gemsec_fb/facebook_clean_data/%s_edges.csv' %s))


GEMSEC_DEEZER_ARGS = ['RO', 'HR', 'HU']
def load_gemsec_deezer(s):
    assert s in GEMSEC_DEEZER_ARGS
    return nx.read_adjlist(os.path.join(DIRNAME, 'snap/deezer_clean_data/%s_edges.csv' % s))


def load_highschool(minutes, scale=1):
    """ Loads the high school dataset with minutes removed 
        (also scales up by networkit's scaling to given scale factor if > 1)
    """
    edges = {}
    with open(os.path.join(DIRNAME, 'snap/highschool/sd02.txt', 'r')) as f:
        for line in f.readlines():
            trip = (line.strip().split('\t'))
            pair = tuple(sorted([int(trip[0]), int(trip[1])]))
            duration = float(trip[-1])
            edges[pair] = edges.get(pair, 0) + duration

    edgelist = [_[0] for _ in edges.items() if _[1] >= (minutes * 3)]
    G = nx.Graph() 
    G.add_edges_from(edgelist)

    if scale == 1:
        return G    
    else:
        raise NotImplementedError("Need to scale up according to LFR Networkit")

ARXIV_COLLAB_ARGS = ['AstroPh', 'CondMat', 'HepPh', 'GrQc', 'HepTh']
def load_arxiv_collabco(s):
    assert s in ARXIV_COLLAB_ARGS
    return nx.read_adjlist(os.path.join(DIRNAME, 'snap/collab/ca-%s.txt' %s)) 





