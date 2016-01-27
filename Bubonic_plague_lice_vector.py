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

class Humans():
    def __init__(self, S, I, R, D):
        "Human parameters"
        self.nu = ((1.04**(1/365.0))-1.0)
        self.mu = (1-(0.96**(1/365.0)))
        self.gamma = 1/26.0
        self.chance_of_recovery = 0.4
        self.S = S
        self.I = I
        self.R = R
        self.D = D
        self.N = self.S + self.I + self.R
        self.SIR = [self.S, self.I, self.R, self.D]
        self.roaming = 0.3
        self.not_roaming = 1.0 - self.roaming
        self.last_exposed = 0
class Lice():
    def __init__(self, S, I, D, bite_rate, transmission_rate_v_to_h, transmission_rate_h_to_v):
        "Lice parameters"
        self.nu = 0.111
        self.gamma = 1/3.0
        self.bite_rate = bite_rate #1.0-5.0
        self.transmission_rate_v_to_h = transmission_rate_v_to_h #0.3
        self.transmission_rate_h_to_v = transmission_rate_h_to_v #1.0
        self.lice_per_human = 10.0
        self.carrying_capacity = None
        self.S = S
        self.I = I
        self.D = D
        self.N = self.S + self.I
        self.SIR = [self.S, self.I, self.D]
        
class Metapop():
    def __init__(self, bite_rate, transmission_rate_v_to_h, transmission_rate_h_to_v)):
        "Metapopulation variables"
        self.neighbors = None
        self.infectious_neighbors_lice = None
        self.N_neighbors = None
        #Parameters to fit
        self.bite_rate = bite_rate
        self.transmission_rate_v_to_h = transmission_rate_v_to_h
        self.transmission_rate_h_to_v = transmission_rate_h_to_v
    def _add_human_lice_objects(self):
        "Healthy human and lice objects"
        self.h = Humans(S = 52.0, I = 0.0, R = 0.0, D = 0.0)
        self.l = Lice(S = 10.0*52.0, I = 0.0, D = 0.0, bite_rate = self.bite_rate, transmission_rate_v_to_h = self.transmission_rate_v_to_h, transmission_rate_h_to_v=self.transmission_rate_h_to_v)
    def _add_infected_humans_and_lice(self):
        "Infected human and lice objects"
        self.h = Humans(S = 51.0, I = 1.0, R = 0.0, D = 0.0)
        self.l = Lice(S = 10.0*51.0, I = 10.0*1.0, D = 0.0, bite_rate = self.bite_rate, transmission_rate_v_to_h = self.transmission_rate_v_to_h, transmission_rate_h_to_v=self.transmission_rate_h_to_v)
    def _human_dynamics_inside_metapop(self):
        "Human equations inside of a cell"
        #Natural deaths
        mu_S = 0
        mu_I = 0
        mu_R = 0

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
        if self.h.S > 0 and self.l.I > 0:
            new_I = min (self.h.S, s.distribution((self.h.not_roaming*self.l.bite_rate*self.l.transmission_rate_v_to_h*(self.l.I/self.h.N)*self.h.S),1))
        if self.h.I > 0:
            new_RD = min (self.h.I, s.distribution((self.h.gamma*self.h.I),1))
        if new_RD > 0:
            new_R = min (new_RD, s.distribution((self.h.chance_of_recovery*new_RD),1))
            new_D = new_RD - new_R

        dS = S_births + R_births - new_I - mu_S
        dI = new_I - new_RD - mu_I
        dR = new_R - mu_R
        dD = new_D + mu_S + mu_I + mu_R

        self.h.S = self.h.S+dS
        self.h.I = self.h.I+dI
        self.h.R = self.h.R+dR
        self.h.D = self.h.D+dD

        self.h.SIR = [self.h.S, self.h.I, self.h.R, self.h.D]
        self.h.N = self.h.S + self.h.I + self.h.R
        self.h.last_exposed = int(new_I)
    def _lice_dynamics_inside_metapop(self):
        "Lice equations inside of a cell"
        self.l.carrying_capacity = self.l.lice_per_human*self.h.N
        new_I = 0
        S_births = 0
        new_D = 0

        if self.l.S > 0 and self.l.carrying_capacity > 0 and (self.l.N/self.l.carrying_capacity) < 1:
            S_births = s.distribution((self.l.nu*self.l.S*(1-(self.l.N/self.l.carrying_capacity))),1)
        if self.h.I > 0 and self.h.N > 0:
            new_I = min (self.l.S, s.distribution((self.l.bite_rate*self.l.transmission_rate_h_to_v*(self.h.I/self.h.N)*self.l.S),1))
        if self.l.I > 0:
            new_D = min (self.l.I, s.distribution((self.l.gamma*self.l.I),1))

        dS = S_births - new_I
        dI = new_I - new_D
        dD = new_D
        
        self.l.S = self.l.S+dS
        self.l.I = self.l.I+dI
        self.l.D = self.l.D+dD

        self.l.SIR = [self.l.S, self.l.I, self.l.D]
        self.l.N = self.l.S + self.l.I
    def _dynamics_between_metapop(self):
        "Equations for humans exposed from neighboring cells"
        new_I = 0
        if self.infectious_neighbors_lice > 0 and self.h.S > 0 and self.N_neighbors:
            new_I = min (self.h.S, s.distribution((self.h.roaming*self.l.bite_rate*self.l.transmission_rate_v_to_h*(self.infectious_neighbors_lice/self.N_neighbors)*self.h.S),1))
        dS = - new_I
        dI = new_I

        self.h.S = self.h.S+dS
        self.h.I = self.h.I+dI

        self.h.SIR = [self.h.S, self.h.I, self.h.R, self.h.D]
        self.h.N = self.h.S + self.h.I + self.h.R
        self.h.last_exposed = self.h.last_exposed + int(new_I)
class Community():
    def __init__(self, end_time, length, width, bite_rate, transmission_rate_v_to_h, transmission_rate_h_to_v):
        "Community variables"
        #General
        self.infected_metapops = 1
        self.time_step = 1.0
        self.end_time = end_time
        self.length = length
        self.width = width
        self.size = self.length*self.width
        self.community = list(None for i in xrange(self.size))
        #Parameters to fit
        self.bite_rate = bite_rate
        self.transmission_rate_v_to_h = transmission_rate_v_to_h
        self.transmission_rate_h_to_v = transmission_rate_h_to_v
        
        #Human
        self.metapop_human_SIR = list(None for i in xrange(self.size))
        self.epidemic_solution_human = None
        self.initial_size = self.size * 52.0
        self.last_exposed_total = 0

        #Lice
        self.metapop_lice_SIR = list(None for i in xrange(self.size))
        self.epidemic_solution_lice = None
        
        #Epidemic
        self.epidemic_duration = 0
        self.distance = None
    def _add_metapops(self):
        "Adds metapop objects to the community variable"
        for i in xrange(len(self.community)):
            self.community[i] = Metapop(bite_rate = self.bite_rate, transmission_rate_v_to_h = self.transmission_rate_v_to_h, transmission_rate_h_to_v=self.transmission_rate_h_to_v)
            self.community[i]._add_human_lice_objects()
    def _add_plague(self):
        "Adds infected metapop objects to the community variable"
        if self.infected_metapops == 0:
            pass
        else:
            for x in xrange(self.infected_metapops):
                i = np.random.randint(len(self.community), size = 1)
                self.community[i]._add_infected_humans_and_lice()
    def _SIR_numbers(self):
        "Adds SIR numbers for each human and lice object to get community totals"
        for i in xrange(len(self.metapop_human_SIR)):
            self.metapop_human_SIR[i] = self.community[i].h.SIR
        self.community_human_totals = list(sum(col) for col in zip(*self.metapop_human_SIR))
        for i in xrange(len(self.metapop_lice_SIR)):
            self.metapop_lice_SIR[i] = self.community[i].l.SIR
        self.community_lice_totals = list(sum(col) for col in zip(*self.metapop_lice_SIR))
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
            self.community[i].infectious_neighbors_lice = sum(self.community[x].l.I for x in self.community[i].neighbors)
            self.community[i].N_neighbors = sum(self.community[x].h.N for x in self.community[i].neighbors)
    def _update_community(self):
        "Updates SIR numbers for calculations"
        map(lambda x:x._human_dynamics_inside_metapop(), self.community), self._SIR_numbers()
        map(lambda x:x._lice_dynamics_inside_metapop(), self.community), self._SIR_numbers()  
        self._calc_neighbors(), map(lambda x:x._dynamics_between_metapop(), self.community), self._SIR_numbers(), self._last_exposure()   
    def epidemic(self):
        "Creates a community of metapopulation objects"
        self._add_metapops(), self._add_plague(), self._SIR_numbers(), self._neighbors()
        self.epidemic_solution_human = [self.community_human_totals]
        self.epidemic_solution_lice = [self.community_lice_totals]

        t = np.linspace(0.0, self.end_time, (self.end_time+1/self.time_step))

        for x in t:
            self._update_community()
            self.epidemic_solution_human.append(self.community_human_totals)
            self.epidemic_solution_lice.append(self.community_lice_totals)
            if self.last_exposed_total > 0:
                self.epidemic_duration = int(x)
            if (x-2*26.0) >= self.epidemic_duration and self.infected_metapops > 0:
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
        data = [int(self.initial_size), int(self.epidemic_duration), self.distance, int(self.deaths/self.initial_size*100), self.bite_rate, self.transmission_rate_v_to_h, self.transmission_rate_h_to_v]
        jsondata = json.dumps(data)
        s.output.write(jsondata)
        if self.distance > 0.0 and int(self.deaths/self.initial_size*100) > 5: 
            s.output.write(jsondata)
        s.outputraw.write(jsondata)
        print data
    def graph(self):
        "Function to call graphing functions"
        self._graph_humans(), self._graph_lice()
    def _graph_humans(self):
        "Exports graph of solution to svg file"
        epidemic_solution_array = np.array(self.epidemic_solution_human)
        human_S_class = epidemic_solution_array[:,0]
        human_I_class = epidemic_solution_array[:,1]
        human_R_class = epidemic_solution_array[:,2]
        human_D_class = epidemic_solution_array[:,3]

        plt.subplot(2,1,1)
        S_line = plt.plot(human_S_class, label='Susceptible', linewidth=2, color='g')
        I_line =plt.plot(human_I_class, label='Infectious', linewidth=2, color='m')
        R_line = plt.plot(human_R_class, label='Recovered', linewidth=2, color='b')
        D_line = plt.plot(human_D_class, label='Dead', linewidth=2, color='r')
        plt.legend(loc=1)
        #plt.xlabel('Days')
        plt.ylabel('Humans')
        plt.title('Bubonic plague with lice without infection')
        plt.grid()
        plt.savefig('lice_model_test.png')  
    def _graph_lice(self):
        epidemic_solution_array = np.array(self.epidemic_solution_lice)
        lice_S_class = epidemic_solution_array[:,0]
        lice_I_class = epidemic_solution_array[:,1]
        lice_D_class = epidemic_solution_array[:,2]

        #Exports graph fo solutution to svg file
        plt.subplot(2,1,2)
        S_line = plt.plot(lice_S_class, label='Susceptible', linewidth=2, color='g')
        I_line =plt.plot(lice_I_class, label='Infectious', linewidth=2, color='m')
        SD_line = plt.plot(lice_D_class, label='Dead', linewidth=2, color='r')
        plt.legend(loc=1)
        plt.xlabel('Days')
        plt.ylabel('Lice')
        plt.grid()
        plt.savefig('lice_model_test.png')
        plt.close()
class Simulator():
    def __init__(self, size_range_min, size_range_max, repeat, randavg):
        self.size_range_min = size_range_min
        self.size_range_max = size_range_max
        self.repeat = repeat
        self.simulation_results = []
        self.randavg = randavg #enable or disable random avg
        self.output = open("licetest1.txt", 'w')
        self.outputraw = open("licetestraw1.txt", 'w')
        self.distance_solution = []
    def distribution(self, lam, size):
      log.debug(' lam %f' % lam)
      log.debug(' size %i' % size)
      if not s.randavg:
        return np.random.poisson(lam, size)[0]
      else:
        return lam
    def simulation(self):
        b=0.3
        th=0.3
        for tv in (0.1, 0.3, 0.5):
            for i in range(self.size_range_min, self.size_range_max+1):
                r=0
                while r < self.repeat:
                    c = Community(width = i, length = i, end_time = 10000.0, bite_rate = b, transmission_rate_v_to_h = tv, transmission_rate_h_to_v=th)
                    c.epidemic()
                    c.output_results()
                if c.distance > 0.0:
                    self.distance_solution.append(c.distance)
                    r=r+1
                    #c.graph()
        s.output.close()
        s.outputraw.close()

if __name__ == "__main__":
    # Parser bits
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

    log = logging.getLogger('Human-Lice')
    log.debug(args.norand)
    #norand = args.norand
    
    size_range_min = 6
    size_range_max = 49
    repeat = 5

    s = Simulator(size_range_min, size_range_max, repeat, args.norand)
    s.simulation()
