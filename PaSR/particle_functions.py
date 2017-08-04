# Particle class to represent the stochatic particles in the PaSR reactor.
# also a lot of the functions for evaluating ... 

import numpy as np
import cantera as ct
import random


class Particle:
    "Stochastic Particles"

    def __init__(self, composition, T, gas):
        
        # set up the composition dict
        self.composition = {}        
        for specie in gas.species():
            if specie.name in composition:
                self.composition[specie.name] = composition[specie.name]
            else:
                self.composition[specie.name] = 0.0

        # set up the Temperature
        self.T = T #[K]
            

    ## Resetting the properties
    def set_temperature(self, T):
        self.T = T

    def set_composition_cantera(self, gas, check_sum=0):
        "Set the composition based on cantera Y list"
        for specie in self.composition:
            # specie here is the string name of a specie
            self.composition[specie] = gas.Y[gas.species_index(specie)]

        # not used frequently but can check if they sum to one
        if check_sum:
            sum = 0.0
            for specie in self.composition:
                sum += self.composition[specie]
            assert (sum == 1.0), "New mass fractions don't sum to one"

    def set_composition(self, new_comp):
        "Set the composition based on a dictionary"
        
        for specie in self.composition:
            if specie in new_comp:
                self.composition[specie] = new_comp[specie]
            else:
                self.composition[specie] = 0.0

            
            

    ## Access property functions
    def get_temperature(self):
        return self.T

    def get_composition(self):
        return self.composition


def inflow_outflow(replacement_particles, oxidizer_inflow_composition,
                   fuel_inflow_composition, oxidizer_inflow_temperature,
                   fuel_inflow_temperature, N_oxidizer_particles):
    "replace the particle with the appropriate inflow compsition and temperature"

    for counter,particle in enumerate(replacement_particles):
        if counter <= N_oxidizer_particles:
            particle.set_composition(oxidizer_inflow_composition)
            particle.set_temperature(oxidizer_inflow_temperature)
        else:
            particle.set_composition(fuel_inflow_composition)
            particle.set_temperature(fuel_inflow_temperature)

def average_properties(particle_list):
    n_particles = len(particle_list) 
    composition_list = [None] * n_particles
    for counter, particle in enumerate(particle_list):
        composition_list[counter] = particle.get_composition()
        
    average_composition = {}
    for specie in composition_list[0]:
        average_composition[specie] = sum(composition[specie] for composition in 
                                          composition_list)/n_particles
        
    average_temperature = 0.0
    for particle in particle_list:
        average_temperature += float(particle.get_temperature()/n_particles)
        
    return average_composition, average_temperature
        
            
def mix_particles(mix_pairs):
    "Takes a N_mix_pairs length list of 2-tuples of particles passed by reference."
    "The particles in each 2-tuple are mixed together with IEM from Pope paper"
    a_list = [random.random() for _ in range(0,len(mix_pairs))]
    for counter,pair in enumerate(mix_pairs):
        a = a_list[counter]
        p1 = mix_pairs[0][0]
        p2 = mix_pairs[0][1]
        ## Mixing compositions
        p1_new_composition = {}
        p2_new_composition = {}
        for specie in p1.get_composition():
            p1_new_composition[specie] = p1.get_composition()[specie] + \
            0.5*a*(p2.get_composition()[specie] -
                   p1.get_composition()[specie])

            p2_new_composition[specie] = p2.get_composition()[specie] + \
            0.5*a*(p1.get_composition()[specie] -
                   p2.get_composition()[specie])

        # now set the particle compositions
        p1.set_composition(p1_new_composition)
        p2.set_composition(p2_new_composition)


        ## Mixing the temperatures
        p1_T = p1.get_temperature()
        p2_T = p2.get_temperature()
        
        p1.set_temperature(p1_T + 0.5*a*(p2_T - p1_T))
        p2.set_temperature(p2_T + 0.5*a*(p1_T - p2_T))

        


