# Particle class to represent the stochatic particles in the PaSR reactor.
# also a lot of the functions for evaluating ...

import numpy as np
import cantera as ct
import random
from scipy.integrate import ode


class Particle:
    "Stochastic Particles"

    def __init__(self, composition, T, gas, Z_oxidizer, Z_fuel,
                 mixture_fraction_element):

        # set up the composition dict
        self.composition = {}
        for specie in gas.species():
            if specie.name in composition:
                self.composition[specie.name] = composition[specie.name]
            else:
                self.composition[specie.name] = 0.0

        # set up the Temperature
        self.T = T #[K]

        # elemental mass fractions of mf_element for inlet streams
        self.Z_oxidizer = Z_oxidizer
        self.Z_fuel = Z_fuel

        # mass fraction element
        self.mixture_fraction_element = mixture_fraction_element

        # initial mixture fraction
        gas.TPY = T, 101000, composition
        self.mixture_fraction = (gas.elemental_mass_fraction(mixture_fraction_element) -
                                 Z_oxidizer)/(Z_fuel - Z_oxidizer)


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

    def get_mixture_fraction(self, gas):
        old_state = gas.TPY
        gas.TPY = self.T, 101000, self.composition
        Z_mf_element = gas.elemental_mass_fraction(self.mixture_fraction_element)
        self.mixture_fraction = (Z_mf_element - self.Z_oxidizer)/ \
                                (self.Z_fuel - self.Z_oxidizer)
        # change the gas back to whatever it was, its passed by reference here
        gas.TPY = old_state
        return self.mixture_fraction

## Functions to investigate the particle list properties

def average_properties(particle_list):
    "returns 2-tuple of arithmetic averages((compostion_dict), temperature)"
    n_particles = len(particle_list)

    # composition
    composition_list = [None] * n_particles
    for counter, particle in enumerate(particle_list):
        composition_list[counter] = particle.get_composition()

    average_composition = {}
    for specie in composition_list[0]:
        average_composition[specie] = sum(composition[specie] for composition in
                                          composition_list)/n_particles
    #temperature
    average_temperature = 0.0
    for particle in particle_list:
        average_temperature += float(particle.get_temperature()/n_particles)

    return average_composition, average_temperature

def mixture_fractions(particle_list, gas):
    "returns a list the same length as particle_list, contains particle mixture fraction"
    mixture_fraction_list = [None] * len(particle_list)
    for counter,particle in enumerate(particle_list):
        mixture_fraction_list[counter] = particle.get_mixture_fraction(gas)
    return mixture_fraction_list

def density(particle_list, gas):
    "returns a list the same length as particle_list, contains particle density"
    density_list = [None] * len(particle_list)
    # remember the state that gas had
    old_state = gas.TPY
    for counter,particle in enumerate(particle_list):
        gas.TPY = particle.get_temperature(), 101000, particle.get_composition()
        _,density_list[counter] = gas.TD
    # set gas back to whatever it was holding earlier
    gas.TPY = old_state
    return density_list

def favre_averaged_mixture_fraction(particle_list, gas):
    density_list = density(particle_list,gas)
    mixture_fraction_list = mixture_fractions(particle_list, gas)
    favre_averaged_mixture_fraction = \
    np.mean([mf*rho for mf,rho in zip(mixture_fraction_list,density_list)]) \
        / np.mean(density_list)
    return favre_averaged_mixture_fraction

def favre_variance_mixture_fraction(particle_list, gas):
    density_list = density(particle_list,gas)
    mixture_fraction_list = mixture_fractions(particle_list, gas)
    famf = favre_averaged_mixture_fraction(particle_list,gas)
    famf_fluctuations_squared = [(mf - famf)**2 for mf in
                         mixture_fraction_list]
    favre_variance_mixture_fraction = \
    np.mean([fluc * rho for fluc,rho in zip(famf_fluctuations_squared,density_list)]) \
    / np.mean(density_list)
    return favre_variance_mixture_fraction
    

## Functions to modify the particle list (mixing and inflow/outflow)

def mix_particles(mix_pairs):
    "Takes a N_mix_pairs length list of 2-tuples of particles passed by reference."
    "The particles in each 2-tuple are mixed together with IEM from Pope paper"
    for counter,pair in enumerate(mix_pairs):
        a = 1.#np.random.uniform()
        p1 = pair[0]
        p2 = pair[1]
        if len(p1.get_composition()) != len(p2.get_composition()):
            print("PROBLEM IN MIXING")
            exit()
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

def mixing_IEM(particle_list,gas, tau_mix, dt):
    #famf = favre_averaged_mixture_fraction(particle_list,gas)

    mean_compositions, T_avg = average_properties(particle_list)
    gas.TPY = gas.T, 101000, mean_compositions
    Y_mean = np.append(gas.Y,T_avg)
    for particle in particle_list:
        gas.TPY = particle.get_temperature(), 101000, particle.get_composition()
        Y = np.append(gas.Y,particle.get_temperature())
        def f(t, Y, Y_mean):
            return [-(Yi - mean)/(2.*tau_mix) for Yi,mean in zip(Y,Y_mean)] 

        r = ode(f).set_integrator('vode', method='adams')
        r.set_initial_value(Y, 0.).set_f_params(Y_mean)
        t1 = dt
        dt_ = dt/50.
        while r.successful() and r.t < t1:
            r.integrate(r.t+dt_)
        new_comp = [y/sum(r.y[0:-1]) for y in r.y[0:-1]]
        new_temp = r.y[-1]
        gas.TPY = new_temp, 101000, new_comp
        particle.set_composition_cantera(gas)
        particle.set_temperature(new_temp)


def inflow_outflow(replacement_particles, oxidizer_inflow_composition,
                   fuel_inflow_composition, oxidizer_inflow_temperature,
                   fuel_inflow_temperature, N_oxidizer_particles):
    "replace the particle_list particles with the appropriate \
    inflow composition and temperature"

    for counter,particle in enumerate(replacement_particles):
        if counter < N_oxidizer_particles:
            particle.set_composition(oxidizer_inflow_composition)
            particle.set_temperature(oxidizer_inflow_temperature)
        else:
            particle.set_composition(fuel_inflow_composition)
            particle.set_temperature(fuel_inflow_temperature)
