# Main PaSR file

import numpy as np
import cantera as ct
import particle_functions as pf
import random
import pickle
import time
import matplotlib.pyplot as plt

### Parameters ###
mech = 'DRM_22.cti' #'r_creck_52.cti' #chemical mechanism
T = 300.0 #[K]
n_particles = 1000

# mass fractions and make sure that the string used here is an actual species
# in the mechanism (same spelling)
initial_composition = {'H2': 0.05, 'O2': 0.2185, 'N2':0.7315}

# time scales
tau_res =  2e-3
tau_mix = 2e-4
dt = 0.1 * min([tau_res, tau_mix]) #[s] should be min(tau_res, tau_mix)/10 
# according to pope
#make sure that dt % cantera_dt = 0
cantera_dt = 1e-5 #[s]

# ratio of mass flow rates Oxidizer/(Fuel + Oxidizer)
P = 0.7

# inlet conditions
oxidizer_inflow_composition = {'O2': 0.23, 'N2':0.77}
fuel_inflow_composition = {'H2' : 0.5, 'N2' : 0.5}
oxidizer_inflow_temperature = 300 #[K]
fuel_inflow_temperature =  300 #[K]

# mixture fraction
mixture_fraction_element = 'N'

# Reaction on/off
reactions = 1

# Restart from previous particle list
restart = 0

## End of settings
#..........................................................................#

# get a mixture object going
gas = ct.Solution(mech)

# equilibrate the intial compposition
gas.TPY = T, 101000, initial_composition
gas.equilibrate('HP')
for counter,mass_fraction in enumerate(gas.Y):
    initial_composition[gas.species_name(counter)] = mass_fraction
T = gas.T

# Mixture fraction setup
gas.TPY = T, 101000, oxidizer_inflow_composition
#elemental mass fraction of mt_element in oxidizer stream
Z_oxidizer = gas.elemental_mass_fraction(mixture_fraction_element)
_, rho_fuel = gas.TD

gas.TPY = T, 101000, fuel_inflow_composition
#elemental mass fraction of mt_element in oxidizer stream
Z_fuel = gas.elemental_mass_fraction(mixture_fraction_element)
_, rho_oxidizer = gas.TD

# create list to hold the particles
if restart:
    particle_list = pickle.load(open('saved_particles.p', 'rb'))
else:
    particle_list = [None] * n_particles
    for i in range(n_particles):
        particle_list[i] = pf.Particle(initial_composition, T, gas, Z_oxidizer,
                                       Z_fuel, mixture_fraction_element)
# convergence switch
converged = 0
old_temperature = pf.average_properties(particle_list)[1]
old_old_temperature = old_temperature
old_composition = pf.average_properties(particle_list)[0]

# total PaSR simulation time
total_simulation_time = 0.


## MAIN SIMULATION LOOP
while (not converged):

    ## Step 1: Inflow/Outflow
    N_replace = round(n_particles * float(dt/tau_res))
    print("Number of pairs for inflow replacement: {}".format(N_replace))
    #how many of them are replaced with oxidizer particle or fuel particle
    #(probably could do better here by spliting the last particle ...
    #rather than rounding the N_replace)
    N_oxidizer_particles = round(P * N_replace)

    replacement_particles = [None] * N_replace
    for counter, _ in enumerate(replacement_particles):
        # get by reference the particles to be replaced with new inflow particles
        # this is random selection with replacement per Pope
        replacement_particles[counter] = random.choice(particle_list)

    pf.inflow_outflow(replacement_particles, oxidizer_inflow_composition,
                      fuel_inflow_composition, oxidizer_inflow_temperature,
                      fuel_inflow_temperature, N_oxidizer_particles)

    ## Step 2: Mixing
    # N_mix_pairs = round(n_particles * float(dt/tau_mix))
    # print("Number of pairs for mixing: {}".format(N_mix_pairs))
    # mix_pairs = [None] * N_mix_pairs
    # for counter,_ in enumerate(mix_pairs):
    #     # randomly select two particles with replacement per Pope
    #     p1 = random.choice(particle_list)
    #     p2 = random.choice(particle_list)
    #     mix_pairs[counter] = (p1,p2)
    # pf.mix_particles(mix_pairs)

    pf.mixing_IEM(particle_list, gas, tau_mix, dt)


    ## Step 3: Reaction
    if reactions:
        for counter,particle in enumerate(particle_list):
            # set gas object to the state of this particle
            gas.TPY = particle.get_temperature(),101000,particle.get_composition()

            # stupidly need to reset the reactor/reactornet to make this work
            r = ct.IdealGasConstPressureReactor(gas)
            sim = ct.ReactorNet([r])

            # particle reaction time i.e. cantera reaction time
            r_time = 0.0
            for _ in range(round(dt/cantera_dt)):
                r_time += cantera_dt
                sim.advance(r_time)

            # update the particles after reaction
            particle.set_composition_cantera(gas)
            particle.set_temperature(r.thermo.T)


    ## Check for convergence

    #composition threshold
    new_composition = pf.average_properties(particle_list)[0]
    composition_difference =sum([abs(new_composition[specie] - old_composition[specie])
                              for specie in new_composition])
    print ("Composition difference: {:2.8f}".format(composition_difference))
    old_composition = new_composition

    #temperature threshold
    temperature_difference1 = abs(pf.average_properties(particle_list)[1]
                                  - old_temperature)
    temperature_difference2 = abs(pf.average_properties(particle_list)[1]
                                  - old_old_temperature)

    if temperature_difference1 < 0.1 \
       and temperature_difference2 < 0.1 \
       and total_simulation_time > 0.1 \
       and composition_difference < 1e-6 \
       and temperature_difference1 + temperature_difference2 + \
       composition_difference > 0.0:
        print("Simultion converged, T difference: {}; {}".format(
            temperature_difference1, temperature_difference2))
        print("Composition difference: {:1.8f}".format(composition_difference))
        converged = 1

    else:
        print("Simultion not converged, difference: {}; {}".format(
            temperature_difference1,
            temperature_difference2))
        total_simulation_time += dt
        print("Simulation time: {:06.4f}\n".format(total_simulation_time))
        old_old_temperature = old_temperature
        old_temperature = pf.average_properties(particle_list)[1]

    print ("temperature: {:6.0f}".format(pf.average_properties(particle_list)[1]))
    print ("FMFL {}".format(pf.favre_averaged_mixture_fraction(particle_list,gas)))
    print ("FAVAR: {}".format(pf.favre_variance_mixture_fraction(particle_list,gas)))
    mf2 = pf.mixture_fractions(particle_list,gas)
    print("MF: {}".format(np.mean(mf2)))
    print("Var: {}".format(np.var(mf2)))


# write out the list of particles to a pickle file
pickle.dump(particle_list, open('saved_particles.p', 'wb'))

print (pf.average_properties(particle_list))

print (np.mean(pf.mixture_fractions(particle_list,gas)))

print (pf.favre_averaged_mixture_fraction(particle_list,gas))
