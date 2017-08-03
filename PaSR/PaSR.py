# Main PaSR file

import numpy as np
import cantera as ct
import Particle as pp
import random

### Parameters ###
mech = 'r_creck_52.cti' #chemical mechanism
T = 300 #[K]
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
# p1 = pp.Particle(initial_composition, T, gas)
# print(p1.get_composition())
# p1.set_composition(oxidizer_inflow_composition)
# print(p1.get_composition())
# exit()




# create list to hold the particles
particle_list = [None] * n_particles 
for i in range(n_particles):
    particle_list[i] = pp.Particle(initial_composition, T, gas)


#Step 1: Inflow/Outflow
N_replace = round(n_particles * (dt/tau_res))
N_oxidizer_particles = round(P * N_replace)
N_fuel_particles = N_replace = N_oxidizer_particles
replacement_particles = [None] * N_replace
for counter, _ in enumerate(replacement_particles):
    replacement_particles[counter] = particle_list[random.randint(0,n_particles -1)]

pp.inflow_outflow(replacement_particles, oxidizer_inflow_composition, 
                  fuel_inflow_composition, oxidizer_inflow_temperature,
                  fuel_inflow_temperature, N_oxidizer_particles)



#for particle in particle_list:
    #print(particle.get_composition()['H2O'])

avg_comp = pp.average_properties(particle_list)
print (avg_comp)
exit()


#Step 3: Reaction
for particle in particle_list[0:1]:
    gas.Y = particle.get_composition()
    r = ct.IdealGasConstPressureReactor(gas)
    sim = ct.ReactorNet([r])
    
    time = 0 # local particle reaction time
    for i in range(round(dt/cantera_dt)):
        time += cantera_dt
        sim.advance(time)
    particle.set_composition_cantera(gas)




