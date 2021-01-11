import gym
import math
import random
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import namedtuple, defaultdict
from itertools import count
from PIL import Image
import heapq


import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import torchvision.transforms as T


import EoNlocal.simulation as sim 
import quarantines as Q 

""" 
General setup to figure out when to quarantine using DQN:

Environment: Graph with {S,I,R} on each node 
State/features: Neural network will take in features by a degree-wise approximation
                i.e., inputs will be a vector of length 3N
                [S_1, I_1, R_1, ..., S_N, I_N, R_N] 
                (or maybe some bucketing thereof) 
Actions : {perform quarantine, do nothing} 

Environment: 
    step(action):
        do-nothing: Some I node is selected and either infects a neighbor 
                    or spontaneously recovers 
        perform-quarantine: All I nodes -> R nodes 
                            Some small number of randomly chosen 
                            S nodes get infected 
    terminal state:
        No infected nodes remaining
"""


class EnvironmentSIR:
    def __init__(self, G, tau, gamma, rho, quarantine_marginal=None):
        """ Maintains the state of an SIR model on a networkx graph.

            Basic strategy is to emulate the nonMarkov_SIR model of EoN 
            but allow for discrete steps (should only be 2*len(G) of them at max)
        ARGS: 
            G : networkx graph object we run the epidemic on 
            tau : float - spreading parameter for the SIR model
            gamma : float - recovery parameter for the SIR model 
            rho : float - fraction of nodes who get infected upon restart 
            quarantine_marginal: int -> float - marginal cost of doing the i^th
                                 quarantine (starting with 1). Defaults to 
                                 a 10x^2 cost
        """
        # Record the args 
        self.G = G 
        self.tau = tau 
        self.gamma = gamma 
        self.rho = rho 
        self.quarantine_marginal = quarantine_marginal or (lambda i: -20 * i)


        # And setup the state to use and initialize the infection
        self.reset()


    def reset(self):
        """ Completely resets the environment """

        # vars needed for RL 
        self.max_degree = max(self.G.degree(), key=lambda p: p[1])[1]        
        self.dwise = torch.zeros((3, self.max_degree + 1))
        for _, deg in self.G.degree():
            self.dwise[0][deg] +=1
        self.num_quarantines = 0
        self.quarantine_history = []
        self.reverse_status = {'S': set(self.G.nodes()), 
                               'I': set(), 
                               'R': set()}        

        # vars needed for EoN
        self.status = defaultdict(lambda: 'S')
        self.rec_time= defaultdict(lambda: -1)
        self.pred_inf_time = defaultdict(lambda: float('inf'))      
        self.times = [0]
        self.S = [self.G.order()]
        self.I = [0]
        self.R = [0]
        self.transmissions = []
        self.trans_rate_fxn = lambda x, y: self.tau
        self.rec_rate_fxn = lambda x: self.gamma        
        self.Q = sim.myQueue(float('inf'))

        self._initialize_infection()



    def degreewise_state(self): 
        # state-getter function
        return self.dwise, self.num_quarantines


    def _initialize_infection(self):
        """ Empties the queue and reinfects some nodes. 
            This method is only called after a reset OR after a quarantine

            Returns: (-#new nodes infected, done)
            where done is a boolean that is True if the epidemic is over 
        """

        # Reset the queue's state
        self.Q = sim.myQueue(float('inf'))

        # If no one left to infect, return terminal state
        if len(self.reverse_status['S']) == 0:
            return 0, True

        # Pick a set of initially infected nodes 
        num_infected = max([int(self.rho * len(self.reverse_status['S'])), 1])
        initial_infecteds = random.sample(self.reverse_status['S'], num_infected)
        # ... and for each initially infected node
        start_t = self.times[-1]
        trans_and_rec_time_fxn = sim._trans_and_rec_time_Markovian_const_trans_
        trans_and_rec_time_args=(self.tau, self.rec_rate_fxn)
        for u in initial_infecteds:
            # run the infect method to add future events to the queue
            self.pred_inf_time[u] = start_t
            sim._process_trans_SIR_(start_t, self.G, None, u, self.times, 
                                    self.S, self.I, self.R, self.Q, self.status,
                                    self.rec_time, self.pred_inf_time, 
                                    self.transmissions, 
                                    trans_and_rec_time_fxn, 
                                    trans_and_rec_time_args)
            # And modify the degree-wise state (doing naively/slowly)
            deg_u = self.G.degree(u)
            self.dwise[0][deg_u] -= 1 
            self.dwise[1][deg_u] += 1 

        # And modify the RL state: 
        self.reverse_status['S'] = self.reverse_status['S'].difference(initial_infecteds)
        self.reverse_status['I'] = self.reverse_status['I'].union(initial_infecteds)


        return -len(initial_infecteds), False

    def step(self, action): #DONE
        """ Takes a timestep.
        ARGS:
            action : int -
                    + if 0, proceeds with a timestep as expected for SIR model 
                    + if 1, does a quarantine 
        RETURNS: 
            (reward, done)
            reward is a float of how much reward we got this step 
            done is a boolean: True if no one is infected after this step 
        """ 

        assert action in [0,1] 
        if action == 0:
            
            return self.do_epidemic_step() 
        else:
            return self.do_quarantine_step()

    def run_to_completion(self):
        # Just debugging steps 
        done = False
        accum = 0
        while not done:
            next_accum, done = self.step(0)
            accum += next_accum 
        return accum


    def plot_IR(self):
        tup = Q.TupleSIR(self.times, self.S, self.I, self.R, tmax=float('inf'))
        tup.plot_single()


    def do_epidemic_step(self): #DONE?
        """ Does a step in the SIR model: pops one element from the queue 
            and computes the reward function from that
        """

        #Note that when finally infected, pred_inf_time is correct
        #and rec_time is correct.  
        #So if return_full_data is true, these are correct
        if len(self.Q) == 0: #????
            return 0, True

        t, counter, fxn, args = heapq.heappop(self.Q._Q_)
        reward = 0
        if fxn == sim._process_rec_SIR_: # RECOVERY EVENT
            node = args[0]
            # modify the degreewise state and compute reward            
            deg = self.G.degree[node]
            self.dwise[1][deg] -= 1 
            self.dwise[2][deg] += 1
            self.reverse_status['I'].remove(node)
            self.reverse_status['R'].add(node)

        else: # INFECTION EVENT
            # modify the degreewise state
            node = args[2]
            if self.status[node] == 'S':            
                deg = self.G.degree[node]
                self.dwise[0][deg] -= 1
                self.dwise[1][deg] += 1
                self.reverse_status['S'].remove(node)
                self.reverse_status['I'].add(node)
                reward = -1            
        fxn(t, *args)
        done = (len(self.Q) == 0)
        return reward, done


    def do_quarantine_step(self):
        """ Modify state from performing a quarantine. 
            Need to: 
            - update the state to reflect quarantine 
            - reset Q 
            - compute reward 
            - reinitialize infection 
        """
        self.num_quarantines += 1
        quarantine_reward = self.quarantine_marginal(self.num_quarantines)

        # Loop through all nodes and modify everything
        qtime = self.times[-1]
        self.quarantine_history.append(qtime)
        for u in self.reverse_status['I']:
            # modify the dwise state
            deg_u = self.G.degree[u]
            self.dwise[1][deg_u] -= 1 
            self.dwise[2][deg_u] += 1 

            # modify all EoN state vars 
            sim._process_rec_SIR_(qtime, u, self.times, self.S, self.I, 
                                  self.R, self.status)

        # modify reverse status
        self.reverse_status['R'] |= self.reverse_status['I']
        self.reverse_status['I'] = set() 


        # Reward comes from adding another quarantine 
        # and number of newly infected nodes 
        reinit_reward, done = self._initialize_infection()
        print
        return reinit_reward + quarantine_reward, done
