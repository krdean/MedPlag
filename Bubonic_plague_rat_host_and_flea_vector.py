#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as plt
import json
import pdb
import argparse
import logging


class Rats():
    def __init__(self, S, I, R, D, roaming):
        self.nu = ((4.57**(1/365.0))-1.0)
        self.mu = 0.00068
        self.carrying_capacity = 52.0 
        self.beta = .0641
        self.gamma = 1/5.15
        self.inherited_resistance = 0.5
        self.chance_of_recovery = 0.10
        self.roaming = roaming
        self.not_roaming = 1.0 - self.roaming

        self.S = S
        self.I = I
        self.R = R
        self.D = D
        self.N = self.S + self.I + self.R
        self.SIR = [self.S, self.I, self.R, self.D]
class Fleas():
    def __init__(self, H, S, I, carrying_capacity_per_rat):
        self.nu = ((21**(1/365.0))-1.0)
        self.mu = 1/5.0 #Off host death rate
        self.carrying_capacity = None
        self.carrying_capacity_per_rat = carrying_capacity_per_rat
        self.searching_efficiency = 0.038

        self.H = H #Average number of fleas on host
        self.I = I #Free susceptible fleas
        self.SIR = [self.H, self.I]
class Humans():
    def __init__(self, S, I, R, D):
        self.nu = ((1.04**(1/365.0))-1.0)
        self.mu = (1.0-(0.96**(1/365.0)))
        self.beta = 0.0641
        self.gamma = 1/26.0
        self.chance_of_recovery = 0.33
        self.last_exposed = 0

        self.S = S
        self.I = I
        self.R = R
        self.D = D
        self.N = self.S + self.I + self.R
        self.SIR = [self.S, self.I, self.R, self.D]
class Metapop():
    def __init__(self, roaming, carrying_capacity_per_rat):
        "Metapopulation variables"
        self.infectious_rat_deaths = None
        self.neighbors = None
        self.force_to_humans = 0
        self.force_to_rats = 0
        self.neighbor_force_to_rats = 0
        self.roaming = roaming
        self.carrying_capacity_per_rat = carrying_capacity_per_rat
    def _add_host_vector_objects(self):
        "Defines healthy rat, flea, human objects"
        self.r = Rats(S = 52.0, I = 0.0, R = 0.0, D = 0.0, roaming = self.roaming)
        self.f = Fleas(H = self.carrying_capacity_per_rat, S = 0.0, I = 0.0, carrying_capacity_per_rat = self.carrying_capacity_per_rat)
        self.h = Humans(S = 52.0, I = 0.0, R = 0.0, D = 0.0)
    def _add_infected_host_vector_objects(self):
        "Defines infected rat, flea, human objects"
        self.r = Rats(S = 51.0, I =1.0, R = 0.0, D = 0.0, roaming = self.roaming)
        self.f = Fleas(H = self.carrying_capacity_per_rat, S = 0.0, I = self.carrying_capacity_per_rat, carrying_capacity_per_rat = self.carrying_capacity_per_rat)
        self.h = Humans(S = 52.0, I = 0.0, R = 0.0, D = 0.0) 
    def _rat_equations_inside(self):
        "Rat dyniamics inside of a cell"
        mu_S = 0
        mu_I = 0
        mu_R = 0
        
        #Natural deaths
        if self.r.S > 0:
            mu_S = min (self.r.S, s.distribution((self.r.mu*self.r.S),1))
        if self.r.I > 0:
            mu_I = min (self.r.I, s.distribution((self.r.mu*self.r.I),1))
        if self.r.R > 0:
            mu_R = min (self.r.R, s.distribution((self.r.mu*self.r.R),1))
        
        S_births = 0
        nu_R = 0
        nu_R_p = 0
        born_R = 0
        S_R_births = 0 
        new_I = 0
        new_RD = 0
        new_R = 0
        new_D = 0
        self.infectious_rat_deaths = 0
        self.noninfectious_rat_deaths = 0
            
        if self.r.S > 0 and (self.r.N/self.r.carrying_capacity) < 1:
            S_births = s.distribution((self.r.nu*self.r.S*(1-(self.r.N/self.r.carrying_capacity))),1)
        if self.r.R > 0 and (self.r.N/self.r.carrying_capacity) < 1:
            nu_R = s.distribution((self.r.nu*self.r.R*(1-(self.r.N/self.r.carrying_capacity))),1)
            if nu_R > 0:
                born_R = min(nu_R, s.distribution((nu_R*self.r.inherited_resistance),1))
            S_R_births = nu_R - born_R
        if self.r.S > 0 and self.force_to_rats > 0:
            new_I = min (self.r.S, s.distribution((self.r.not_roaming*self.r.beta*(self.r.S/self.r.N)*self.force_to_rats),1))
        if self.r.I > 0:
            new_RD = min (self.r.I, s.distribution((self.r.gamma*self.r.I),1))
            new_R = s.distribution((self.r.chance_of_recovery*new_RD),1)
            new_D = new_RD - new_R
        
        dS = S_births + S_R_births - new_I - mu_S
        dI = new_I - new_RD - mu_I
        dR = born_R + new_R - mu_R
        dD = new_D + mu_S + mu_I + mu_R

        if (new_D + mu_I) > 0:
            self.infectious_rat_deaths = new_D + mu_I
        
        self.r.S = self.r.S + dS
        self.r.I = self.r.I + dI
        self.r.R = self.r.R + dR
        self.r.D = self.r.D + dD
        
        self.r.SIR = [self.r.S, self.r.I, self.r.R, self.r.D]
        self.r.N = self.r.S + self.r.I + self.r.R
        self.f.carrying_capacity = self.r.N*self.f.carrying_capacity_per_rat
    def _flea_equations_inside(self):
        "Flea dynamics inside of a cell"
        flea_growth = 0
        new_I = 0
        mu_I = 0
        self.force_to_humans = 0
        self.force_to_rats = 0
        
        #Growth and decay
        if self.f.H > 0 and self.f.carrying_capacity > 0 and ((self.f.H*self.r.N)/self.f.carrying_capacity) < 1:
            flea_growth = s.distribution((self.f.nu*self.f.H*self.r.N*(1-(self.f.H*self.r.N)/self.f.carrying_capacity)),1)
        elif self.f.H > 0 and self.f.carrying_capacity > 0 and ((self.f.H*self.r.N)/self.f.carrying_capacity) > 1:
            flea_growth = -s.distribution((self.f.nu*self.f.H*self.r.N*abs(1-(self.f.H*self.r.N/self.f.carrying_capacity))),1)
        #Starvation
        if self.f.I > 0:
            mu_I = min (self.f.I, s.distribution((self.f.mu*self.f.I),1)) #Free infectious fleas dying of starvation
        #New free fleas
        if self.infectious_rat_deaths > 0 and self.f.H > 0:
            new_I = int(self.infectious_rat_deaths*self.f.H)
        #Force of infection to humans
        if self.f.I > 0 and self.r.N > 0:
            self.force_to_humans = min(self.f.I, s.distribution(self.f.I*np.exp(-self.f.searching_efficiency*self.r.N),1))
        #Force of infection to rats
            self.force_to_rats = self.f.I-self.force_to_humans
        
        total_H_changes = flea_growth + self.force_to_rats
        if self.r.N > 0:
            avg_H_changes = total_H_changes/self.r.N
        else:
            avg_H_changes = 0

        dH = avg_H_changes
        dI = new_I - mu_I

        self.f.H = self.f.H + dH #Average number of fleas on rats
        self.f.I = self.f.I + dI #Free infectious fleas
        self.f.SIR = [self.f.H, self.f.I]
    def _human_equations_inside(self):
        "Human dynamics inside of a cell"       
        mu_S = 0
        mu_I = 0
        mu_R = 0
        
        #Natural Deaths
        if self.h.S > 0:
            mu_S = min (self.h.S, s.distribution((self.h.mu*self.h.S),1))
        if self.h.I > 0:
            mu_I = min (self.h.I, s.distribution((self.h.mu*self.h.I),1))
        if self.h.R > 0:
            mu_R = min (self.h.R, s.distribution((self.h.mu*self.h.R),1))
        
        S_births = 0
        R_births = 0
        new_I = 0
        new_RD = 0
        new_R = 0
        new_D = 0
        
        if self.h.S > 0:
            S_births = s.distribution((self.h.nu*self.h.S),1)
        if self.h.R > 0:
            R_births = s.distribution((self.h.nu*self.h.R),1)
        if self.h.S > 0 and self.force_to_humans > 0:
            new_I = min (self.h.S, s.distribution((self.h.beta*(self.h.S/self.h.N)*self.force_to_humans),1))
        if self.h.I > 0:
            new_RD = min (self.h.I, s.distribution((self.h.gamma*self.h.I),1))
            new_R = min(new_RD, s.distribution((self.h.chance_of_recovery*new_RD),1))
            new_D = new_RD - new_R 

        dS = S_births + R_births - new_I - mu_S
        dI = new_I - new_RD - mu_I
        dR = new_R - mu_R
        dD = new_D + mu_S + mu_I + mu_R
        
        self.h.S = self.h.S + dS
        self.h.I = self.h.I + dI
        self.h.R = self.h.R + dR
        self.h.D = self.h.D + dD
        
        self.h.SIR = [self.h.S, self.h.I, self.h.R, self.h.D]
        self.h.N = self.h.S + self.h.I + self.h.R
        self.h.last_exposed = int(new_I)
    def _dynamics_between_metapop(self):
        "Equations for rats exposed from neighboring cells"
        new_I = 0
        if self.N_neighbors_rats > 0 and self.neighbor_force_to_rats > 0 and self.r.S > 0:
            new_I = min (self.r.S, s.distribution((self.r.roaming*self.r.beta*(self.r.S/self.N_neighbors_rats)*self.neighbor_force_to_rats),1))
        dS = - new_I
        dI = new_I

        self.r.S = self.r.S+dS
        self.r.I = self.r.I+dI

        self.r.SIR = [self.r.S, self.r.I, self.r.R, self.r.D]
        self.r.N = self.r.S + self.r.I + self.r.R      
        self.f.carrying_capacity = self.r.N*self.f.carrying_capacity_per_rat
class Community():
    def __init__(self, end_time, length, width, roaming, carrying_capacity_per_rat):
        "Community variables"
        #General
        self.infected_metapops = 1
        self.time_step = 1.0
        self.end_time = end_time
        self.length = length
        self.width = width
        self.size = self.length*self.width
        self.community = list(None for i in xrange(self.size))
        self.roaming = roaming
        self.carrying_capacity_per_rat = carrying_capacity_per_rat

        #Human
        self.metapop_human_SIR = list(None for i in xrange(self.size))
        self.epidemic_solution_human = None
        self.initial_size = self.size * 52.0
        self.last_exposed_total = 0

        #Rat
        self.metapop_rat_SIR = list(None for i in xrange(self.size))
        self.epidemic_solution_rats = None

        #Flea
        self.metapop_flea_SIR = list(None for i in xrange(self.size))
        self.epidemic_solution_flea = None
        
        #Epidemic
        self.epidemic_duration = 0
        self.distance = None
        self.human_epidemic_duration = 0     
    def _add_metapops(self):
        "Adds metapop objects to the community variable"
        for i in xrange(len(self.community)):
            self.community[i] = Metapop(roaming = self.roaming, carrying_capacity_per_rat = self.carrying_capacity_per_rat)
            self.community[i]._add_host_vector_objects()
    def _add_plague(self):
        "Adds infected metapop obects to the community variable"
        if self.infected_metapops == 0:
            pass
        else:
            for x in xrange(self.infected_metapops):
                i = np.random.randint(len(self.community), size = 1)
                self.community[i]._add_infected_host_vector_objects()
    def _SIR_numbers(self):
        "Adds SIR numbers for vector and host objects to get community totals"
        for i in xrange(len(self.metapop_rat_SIR)):
            self.metapop_rat_SIR[i] = self.community[i].r.SIR
        self.community_rat_totals = list(sum(col) for col in zip(*self.metapop_rat_SIR))
        for i in xrange(len(self.metapop_flea_SIR)):
            self.metapop_flea_SIR[i] = self.community[i].f.SIR
        self.community_flea_totals = list(sum(col) for col in zip(*self.metapop_flea_SIR))
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
        for i in xrange(len(self.community)):
            self.community[i].neighbor_force_to_rats = sum(self.community[x].force_to_rats for x in self.community[i].neighbors)
            self.community[i].N_neighbors_rats = sum(self.community[x].r.N for x in self.community[i].neighbors)
    def _update_community(self):
        "Updates SIR numbers for calculations"
        map(lambda x:x._rat_equations_inside(), self.community), self._SIR_numbers()
        map(lambda x:x._flea_equations_inside(), self.community), self._SIR_numbers()
        map(lambda x:x._human_equations_inside(), self.community), self._SIR_numbers()
        self._calc_neighbors(), map(lambda x:x._dynamics_between_metapop(), self.community), self._SIR_numbers(), self._last_exposure()
    def epidemic(self):
        "Creates a community of metapopulation objects"
        self._add_metapops(), self._add_plague(), self._SIR_numbers(), self._neighbors()
        self.epidemic_solution_rat = [self.community_rat_totals]
        self.epidemic_solution_flea = [self.community_flea_totals]
        self.epidemic_solution_human = [self.community_human_totals]
        
        t = np.linspace(0.0, self.end_time, (self.end_time+1/self.time_step))
        for x in t:
            self._update_community()
            self.epidemic_solution_rat.append(self.community_rat_totals)
            self.epidemic_solution_flea.append(self.community_flea_totals)
            self.epidemic_solution_human.append(self.community_human_totals)
            if self.last_exposed_total > 0:
                self.epidemic_duration = int(x)
            if (x-2*26.0) >= self.epidemic_duration and self.infected_metapops > 0: #26 avg latent + infectious period
                break
        human_solution_array = np.array(self.epidemic_solution_human)
        I_human_list = human_solution_array[:,1]
        first_exposed = np.argmax(I_human_list>0)
        if self.epidemic_duration > 0:
            self.human_epidemic_duration = self.epidemic_duration - (first_exposed-1)
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
        self.pointDist(self.initial_size, self.human_epidemic_duration)
        data = [int(self.initial_size), int(self.epidemic_duration), int(self.human_epidemic_duration), self.distance, int(self.deaths), int(self.deaths/self.initial_size*100), self.roaming, self.carrying_capacity_per_rat]
        jsondata = json.dumps(data)
        if self.distance > 0.0: 
            s.output.write(jsondata)
        s.outputraw.write(jsondata)
        print data
    def graph(self):
        self._graph_human()
        self._graph_rat()
        self._graph_flea()
    def _graph_human(self):
        "Exports graph of solution to svg file"
        #Makes array with columns that are made into new arrays for plotting classes
        epidemic_solution_array_human = np.array(self.epidemic_solution_human)
        human_S_class = epidemic_solution_array_human[:,0]
        human_I_class = epidemic_solution_array_human[:,1]
        human_R_class = epidemic_solution_array_human[:,2]
        human_D_class = epidemic_solution_array_human[:,3]

        #Exports graph for solutution to svg file
        plt.subplot(4,1,1)
        S_line = plt.plot(human_S_class, label='Susceptible', linewidth=2, color='g')
        I_line =plt.plot(human_I_class, label='Infectious', linewidth=2, color='m')
        R_line = plt.plot(human_R_class, label='Recovered', linewidth=2, color='b')
        D_line = plt.plot(human_D_class, label='Dead', linewidth=2, color='r')
        plt.legend(loc=1)
        plt.xlabel('Days')
        plt.ylabel('Humans')
        plt.title('Bubonic with Rats and Fleas without infection')
        plt.grid()
        plt.savefig('rat_model_no_plague.svg')
    def _graph_rat(self):
        #Makes array with columns that are made into new arrays for plotting classes
        epidemic_solution_array_rat = np.array(self.epidemic_solution_rat)
        rat_S_class = epidemic_solution_array_rat[:,0]
        rat_I_class = epidemic_solution_array_rat[:,1]
        rat_R_class = epidemic_solution_array_rat[:,2]
        rat_D_class = epidemic_solution_array_rat[:,3]

        #Exports graph for solutution to svg file
        plt.subplot(4,1,2)
        rat_S_line = plt.plot(rat_S_class, label='Susceptible', linewidth=2, color='g')
        rat_I_line =plt.plot(rat_I_class, label='Infectious', linewidth=2, color='m')
        rat_R_line = plt.plot(rat_R_class, label='Recovered', linewidth=2, color='b')
        rat_D_line = plt.plot(rat_D_class, label='Dead', linewidth=2, color='r')
        plt.legend(loc=1)
        plt.xlabel('Days')
        plt.ylabel('Rats')
        plt.grid()
        plt.savefig('rat_model_no_plague.svg')
    def _graph_flea(self):
        epidemic_solution_array_flea = np.array(self.epidemic_solution_flea)
        flea_average = epidemic_solution_array_flea[:,0]
        flea_I_class = epidemic_solution_array_flea[:,1]

        #Exports graph for solutution to svg file
        plt.subplot(4,1,3)
        flea_I_line = plt.plot(flea_I_class, label='Free infectious', linewidth=2, color='y')
        plt.legend(loc=1)
        plt.xlabel('Days')
        plt.ylabel('Fleas')
        plt.grid()
        plt.savefig('rat_model_no_plague.svg')

        plt.subplot(4,1,4)
        flea_average = plt.plot(flea_average/self.size, label='Average per rat', linewidth=2, color='b')
        plt.legend(loc=1)
        plt.xlabel('Days')
        plt.ylabel('Flea index')
        plt.grid()
        plt.savefig('rat_model_no_plague.svg')

class Simulator():
    def __init__(self, size_range_min, size_range_max, repeat, randavg):
        self.size_range_min = size_range_min
        self.size_range_max = size_range_max
        self.repeat = repeat
        self.simulation_results = []
        self.randavg = randavg #enable or disable random avg
        self.output = open("rat_simulation_outputtest.txt", 'w')
        self.outputraw = open("rat_simulation_outputtestraw.txt", 'w')
        self.distance_solution = []  
    def distribution(self, lam, size):
        log.debug(' lam %f' % lam)
        log.debug(' size %i' % size)
        if not s.randavg:
            return np.random.poisson(lam, size)[0]
        else:
            return lam  
    def simulation(self):
        m = 0.3
        K_f =  4.0 
        for i in range(self.size_range_min, self.size_range_max+1):
            r = 0
            while r < self.repeat: 
                c = Community(width = i, length = i, end_time = 10000.0, roaming = m, carrying_capacity_per_rat = K_f)
                c.epidemic()
                c.output_results()
                if c.distance > 0.0:      
                    self.distance_solution.append(c.distance)
                    r = r+1
                c.graph()
        #print np.mean(self.distance_solution)
        self.output.close()  
        s.outputraw.close()        
if __name__ == "__main__":
    # Parser bits
    parser = argparse.ArgumentParser(description='Rat_flea_human Simulator')
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
    
    log = logging.getLogger('Rat_flea_human')
    log.debug(args.norand)
    
    #Run the simulations
    size_range_min = 10
    size_range_max = 10
    repeat = 1
    s = Simulator(size_range_min, size_range_max, repeat, args.norand)
    s.simulation()
