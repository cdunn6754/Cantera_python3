# Main PaSR file

import numpy as np
import cantera as ct
import particle_functions as pf
import random
import pickle

### Parameters ###
mech = 'r_creck_52.cti' #chemical mechanism
T = 950.0 #[K]
n_particles = 1000
# mass fractions and make sure that the string used here is an actual species
# in the mechanism (same spelling)
initial_composition = {'H2': 0.05, 'O2': 0.2185, 'N2':0.7315}
tau_res =  1
tau_mix = 0.05
dt = 1e-2 #[s]
#make sure that dt % cantera_dt = 0
cantera_dt = 1e-4 #[s]
#ratio of mass flow rates Oxidizer/(Fuel + Oxidizer)
P = 0.7
oxidizer_inflow_composition = {'O2': 0.23, 'N2':0.77}
fuel_inflow_composition = {'H2' : 0.5, 'N2' : 0.5}
oxidizer_inflow_temperature = 300 #[K]
fuel_inflow_temperature =  300 #[K]



#..........................................................................#

# get a mixture object going
gas = ct.Solution(mech)

## TESTING ###

# mix testing
# p1 = pf.Particle(initial_composition, 500, gas)
# p2 = p.Particle(oxidizer_inflow_composition, T, gas)
# print(p2.get_composition()['N2'] + p1.get_composition()['N2'])
# pf.mix_particles( [(p1,p2)] )
# print(p2.get_composition()['N2'] + p1.get_composition()['N2'])
# print(p2.get_temperature() + p1.get_temperature())


# exit()




# create list to hold the particles
particle_list = [None] * n_particles 
for i in range(n_particles):
    particle_list[i] = pf.Particle(initial_composition, T, gas)

converged = 0
while (not converged):
    ## Step 1: Inflow/Outflow
    N_replace = round(n_particles * float(dt/tau_res))
    print("Number of pairs for inflow replacement: {}".format(N_replace))
    #how many of them are replaced with oxidizer particle or fuel particle
    N_oxidizer_particles = round(P * N_replace) 
    N_fuel_particles = N_replace = N_oxidizer_particles
    replacement_particles = [None] * N_replace
    for counter, _ in enumerate(replacement_particles):
        # get by reference the particles to be replaced with new inflow particles
        # this is random selection with replacement per Pope
        replacement_particles[counter] = particle_list[random.randint(0,n_particles -1)]

    pf.inflow_outflow(replacement_particles, oxidizer_inflow_composition, 
                      fuel_inflow_composition, oxidizer_inflow_temperature,
                      fuel_inflow_temperature, N_oxidizer_particles)
    ## Step 2: Mixing
    N_mix_pairs = round(n_particles * float(dt/tau_mix))
    print("Number of pairs for mixing: {}".format(N_mix_pairs))
    mix_pairs = [None] * N_mix_pairs
    for counter,_ in enumerate(mix_pairs):
        # randomly select two particles with replacement per pope
        p1 = particle_list[random.randint(0,n_particles - 1)]
        p2 = particle_list[random.randint(0,n_particles - 1)]
        mix_pairs[counter] = (p1,p2)
    pf.mix_particles(mix_pairs)
    
    for particle in particle_list:
        for specie in particle.get_composition():
            if particle.get_composition()[specie] < 0:
                print (particle.get_composition()[specie])

    
    r = ct.IdealGasConstPressureReactor(gas)
    sim = ct.ReactorNet([r])

    ## Step 3: Reaction
    for particle in particle_list:
        #gas.TPY = particle.get_temperature(),101000,particle.get_composition()
        sim.set_initial_time(0.0) # need to reset integrator
        r.thermo.TPY = particle.get_temperature(),101000,particle.get_composition()
        r.syncState()
        sim.reinitialize()
        #print(particle.get_temperature())
        #print(r.thermo.T)
        time = 0.0 # local particle reaction time
        for _ in range(round(dt/cantera_dt)):
            time += cantera_dt
            sim.advance(time)
        particle.set_composition_cantera(gas)
        particle.set_temperature(r.thermo.T)
    converged =1

pickle.dump(particle_list, open('saved_particles.p', 'wb'))
