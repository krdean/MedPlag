#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['font.family']='serif'
import matplotlib.pyplot as plt
import json
import pdb
import argparse
import logging

#pdb.set_trace()
#python -m pdb script.py

class Humans():
    def __init__(self, S, E, I, D, roaming, beta):
        "Human variables"
        self.nu = ((1.04**(1/365.0))-1.0)
        self.mu = (1-(0.96**(1/365.0)))
        self.beta = beta
        self.sigma = 1/4.3
        self.gamma = 1/2.5
        self.S = S
        self.E = E      
        self.I = I
        self.D = D
        self.N = self.S + self.E + self.I
        self.SIR = [self.S, self.E, self.I, self.D]
        self.last_exposed = 0
        self.roaming = roaming
        self.not_roaming = 1.0 - self.roaming
class Metapop():
    def __init__(self, roaming, beta):
        "Metapopulation variables"
        self.neighbors = None
        self.infectious_neighbors = None
        self.N_neighbors = None
        self.roaming = roaming
        self.beta = beta
    def _add_human_objects(self):
        "Defines healthy human objects"
        self.h = Humans(S = 52.0, E = 0.0, I = 0.0, D = 0.0, roaming = self.roaming, beta = self.beta)
    def _add_infected_humans(self):
        "Defines infected human objects"
        self.h = Humans(S = 51.0, E = 1.0, I = 0.0, D = 0.0, roaming = self.roaming, beta = self.beta)
    def _human_dynamics_inside_metapop(self):
        "Calculates and updates S, E, I, D classes"
        
        #Calculates people changing classes due to births and plague
        S_births = 0
        new_E = 0
        new_I = 0
        new_D = 0

        #Births
        if self.h.S > 0:
            S_births = s.distribution((self.h.nu*self.h.S),1)

        #Plague
        if self.h.S > 0 and self.h.I > 0:
            new_E = min (self.h.S, s.distribution((self.h.not_roaming*self.h.beta*self.h.S*self.h.I/self.h.N),1))
        if self.h.E > 0:
            new_I = min (self.h.E, s.distribution((self.h.sigma*self.h.E),1))
        if self.h.I > 0:
            new_D = min (self.h.I, s.distribution((self.h.gamma*self.h.I),1))
        
        #Natural Deaths
        mu_S = 0
        mu_E = 0
        mu_I = 0

        if self.h.S > 0:
            mu_S = min (self.h.S, s.distribution((self.h.mu*self.h.S),1))
        if self.h.E > 0:
            mu_E = min (self.h.E, s.distribution((self.h.mu*self.h.E),1))
        if self.h.I > 0:
            mu_I = min (self.h.I, s.distribution((self.h.mu*self.h.I),1))

        #Sums changes for each class
        dS = S_births - new_E - mu_S
        dE = new_E - new_I - mu_E
        dI = new_I -new_D - mu_I
        dD = new_D + mu_S + mu_E + mu_I

        #Updates changes
        self.h.S = self.h.S + dS
        self.h.E = self.h.E + dE
        self.h.I = self.h.I + dI
        self.h.D = self.h.D + dD
        
        self.h.N = self.h.S + self.h.E + self.h.I
        self.h.SIR = [self.h.S, self.h.E, self.h.I, self.h.D]
        self.h.last_exposed = int(new_E)
    def _dynamics_between_metapop(self):
        "Equations for humans exposed from neighboring cells"
        #Plague exposure from neighbors
        new_E = 0
        if self.infectious_neighbors > 0 and self.h.S > 0 and self.N_neighbors > 0:
            new_E = min (self.h.S, s.distribution((self.h.roaming*self.h.beta*self.h.S*self.infectious_neighbors/self.N_neighbors),1))
        dS = - new_E
        dE = new_E

        #Updates changes
        self.h.S = self.h.S + dS
        self.h.E = self.h.E + dE

        self.h.N = self.h.S + self.h.E + self.h.I
        self.h.SIR = [self.h.S, self.h.E, self.h.I, self.h.D]
        self.h.last_exposed = self.h.last_exposed + int(new_E)
class Community():
    def __init__ (self, end_time, length, width, roaming, beta):
        "Community variables"
        #General
        self.infected_metapops = 1
        self.time_step = 1.0 #days
        self.end_time = end_time
        self.length = length
        self.width = width
        self.size = self.length*self.width
        self.community = list(None for i in xrange(self.size))
        self.roaming = roaming
        self.beta = beta
        
        #Human
        self.metapop_human_SIR = list(None for i in xrange(self.size))
        self.epidemic_solution_human = None
        self.initial_size = self.size * 52.0
        self.last_exposed_total = 0

        #Epidemic
        self.epidemic_duration = 0
        self.distance = None
    def _add_metapops(self):
        "Adds metapop objects to the community variable"
        for i in xrange(len(self.community)):
            self.community[i] = Metapop(roaming = self.roaming, beta = self.beta)
            self.community[i]._add_human_objects()
    def _add_plague(self):
        "Adds infected metapop objects to the community variable"
        if self.infected_metapops == 0:
            pass
        else:
            for x in xrange(self.infected_metapops):
                i = np.random.randint(len(self.community), size = 1)
                self.community[i]._add_infected_humans()
    def _SIR_numbers(self):
        "Adds SIR numbers for each human object to get community totals"
        for i in xrange(len(self.metapop_human_SIR)):
            self.metapop_human_SIR[i] = self.community[i].h.SIR
        self.community_human_totals = list(sum(col) for col in zip(*self.metapop_human_SIR))
    def _last_exposure(self):
        self.last_exposed_total = 0
        for i in xrange(self.size):
            self.last_exposed_total = self.last_exposed_total + self.community[i].h.last_exposed
    def _neighbors(self):
        "Identifies neighboring cells"
        for i in xrange(len(self.community)):
            if i == 0: #TOP LEFT
                self.community[i].neighbors = [i+1, i+self.width]
            elif 1 <= i <= (self.width-2): #FIRST ROW
                self.community[i].neighbors = [i-1, i+1, i+self.width]
            elif i == (self.width-1): #TOP RIGHT
                self.community[i].neighbors = [i-1, i+self.width]
            elif i == ((self.width*self.length)-1): #BOTTOM RIGHT
                self.community[i].neighbors = [i-1, i-self.width]
            elif ((self.width*self.length)-(self.width-1)) <= i <= ((self.width*self.length)-1): #LAST ROW
                self.community[i].neighbors = [i-1, i+1, i-self.width]
            elif i == ((self.width*self.length)-(self.width)): #BOTTOM LEFT
                self.community[i].neighbors = [i+1, i-self.width]
            elif i%(self.width) == 0: #LEFT EDGE
                self.community[i].neighbors = [i+1, i+self.width, i-self.width]
            elif (i-(self.width-1))%(self.width) == 0: #RIGHT EDGE
                self.community[i].neighbors = [i-1, i+self.width, i-self.width]
            else: #CENTER
                self.community[i].neighbors = [i-1, i-self.width, i+1, i+self.width] 
    def _calc_neighbors(self):
        "Calculates the number of infectious people is nearby cells"
        #Defines variables infectious_neighbors and N_neighbors
        for i in xrange(len(self.community)):
            self.community[i].infectious_neighbors = sum(self.community[x].h.I for x in self.community[i].neighbors)
            self.community[i].N_neighbors = sum(self.community[x].h.N for x in self.community[i].neighbors)
    def _update_community(self):
        "Updates SIR numbers for calculations"
        map(lambda x:x._human_dynamics_inside_metapop(), self.community), self._SIR_numbers()
        self._calc_neighbors(), map(lambda x:x._dynamics_between_metapop(), self.community), self._SIR_numbers(), self._last_exposure()
    def epidemic(self):
        "Creates a community of metapopulation objects"
        self._add_metapops(), self._add_plague(), self._SIR_numbers(), self._neighbors()
        self.epidemic_solution_human = [self.community_human_totals]

        t = np.linspace(0.0, self.end_time, (self.end_time+1/self.time_step))
        for x in t:
            self._update_community()
            self.epidemic_solution_human.append(self.community_human_totals)
            if self.last_exposed_total > 0:
                self.epidemic_duration = int(x)
            if (x-2*6.8) >= self.epidemic_duration and self.infected_metapops > 0: #6.8 avg latent + infectious period
                break
    def pointDist(self, x0, y0):
        if y0 > 0:
            x0=x0/1000.0 #duration (months)
            y0=y0/30.0 #pop size (thousands)
            y=3.031+(.132*x0) #olea line
            pointDist=abs(y-y0) #distance
            self.distance = pointDist
        else:
            self.distance = None
    def output_results(self):
        self.deaths = int(self.epidemic_solution_human[self.epidemic_duration][3])
        self.pointDist(self.initial_size, self.epidemic_duration)
        data = [int(self.initial_size), int(self.epidemic_duration), self.distance, int(self.deaths/self.initial_size*100), self.roaming, self.beta]    
        jsondata = json.dumps(data)
        if self.distance > 0.0: 
            s.output.write(jsondata)
        s.outputraw.write(jsondata)
        print data
    def graph_humans(self):
        "Exports graph of solution to png file"
        epidemic_solution_array = np.array(self.epidemic_solution_human)
        self.human_S_class = epidemic_solution_array[:,0]
        self.human_E_class = epidemic_solution_array[:,1]
        self.human_I_class = epidemic_solution_array[:,2]
        self.human_D_class = epidemic_solution_array[:,3]

        #Exports graph fo solutution to png file
        S_line = plt.plot(self.human_S_class, label='Susceptible', linewidth=1, color='g')
        E_line =plt.plot(self.human_E_class, label='Exposed', linewidth=1, color='c')
        I_line = plt.plot(self.human_I_class, label='Infectious', linewidth=1, color='m')
        D_line = plt.plot(self.human_D_class, label='Dead', linewidth=1, color='r')
        plt.subplot(111)
        #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
        plt.xlabel('Days')
        plt.ylabel('People')
        #plt.title('Pneumonic plague model without infection')
        #plt.grid()
        plt.minorticks_on()
        plt.savefig('human_pneumonic_lowbeta.png', bbox_inches='tight')

class Simulator():
    def __init__(self, size_range_min, size_range_max, repeat, randavg):
        self.size_range_min = size_range_min
        self.size_range_max = size_range_max
        self.repeat = repeat
        self.randavg = randavg #enable or disable random avg
        self.output = open("pneumonic_bestfit5.txt", 'w')
        self.outputraw = open("pneumonic_bestfitraw5.txt", 'w')
        self.distance_solution = []
    def distribution(self, lam, size):
        log.debug(' lam %f' % lam)
        log.debug(' size %i' % size)
        if not s.randavg:
            return np.random.poisson(lam, size)[0]
        else:
            return lam
    def simulation(self):
        for b in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
            for m in (0.2, 0.3, 0.4):
                for i in range(self.size_range_min, self.size_range_max+1):
                    r=0
                    while r < self.repeat:
                        c = Community(width = i, length = i, end_time = 10000.0, roaming = m, beta = b)
                        c.epidemic()
                        c.output_results()
                        if c.distance > 0.0 and int(c.deaths/c.initial_size*100) > 5:      
                            self.distance_solution.append(c.distance)
                            r = r+1
                            #c.graph_humans()
        print np.mean(self.distance_solution)
        s.output.close()
        s.outputraw.close()

if __name__ == "__main__":
    #Parser bits
    parser = argparse.ArgumentParser(description='Human Lice Simulator')
    parser.add_argument('--debug', dest='debug', action='store_true',
                        default=False, help='enable debug mode')
    parser.add_argument('--norandom', dest='norand', action='store_true',
                        default=False, help='disables randomness')
    args = parser.parse_args()

    #Build logging
    if args.debug:
      logging.basicConfig(level=logging.DEBUG)
    else:
      logging.basicConfig(level=logging.INFO)

    log = logging.getLogger('Human-pneumonic')
    log.debug(args.norand)
    #norand = args.norand

    size_range_min = 6
    size_range_max = 49
    repeat = 5
    
    s = Simulator(size_range_min, size_range_max, repeat, args.norand)
    s.simulation()

