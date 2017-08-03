# Particle class to represent the stochatic particles in the PaSR reactor.

import numpy as np
import cantera as ct


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
            

    # Resetting the properties
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

    def set_composition(self, new_comp, check_sum=0):
        "Set the composition based on a dictionary"
        for specie in self.composition:
            self.composition[specie] = 0.0
            if specie in new_comp:
                self.composition[specie] = new_comp[specie]
                
            
            

    # Access functions
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
        
            
