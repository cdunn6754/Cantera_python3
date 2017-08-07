import numpy as np
import cantera as ct
import particle_functions as pf
import random
import pickle
import time
#import matplotlib.pyplot as plt

### Parameters ###
mech = 'DRM_22.cti' #'r_creck_52.cti' #chemical mechanism
T = 950.0 #[K]
n_particles = 1000

# mass fractions and make sure that the string used here is an actual species
# in the mechanism (same spelling)
initial_composition = {'H2': 0.05, 'O2': 0.2185, 'N2':0.7315}

# time scales
tau_res =  2e-3
tau_mix = 1e-3
#[s] should be min(tau_res, tau_mix)/10 according to pope
dt = 0.1 * min([tau_res,tau_mix]) 
#make sure that dt % cantera_dt = 0
cantera_dt = 1e-5 #[s]

# ratio of mass flow rates Oxidizer/(Fuel + Oxidizer)
P = 0.7

# inlet conditions
oxidizer_inflow_composition = {'O2': 0.23, 'N2':0.77}
fuel_inflow_composition = {'H2' : 0.5, 'N2' : 0.5}
oxidizer_inflow_temperature = 300 #[K]
fuel_inflow_temperature =  300 #[K]
inflow_key_names = ['O2', 'N2', 'H2']

# get a mixture object going
gas = ct.Solution(mech)

# mixture fraction
mixture_fraction_element = 'O'

# Reaction on/off
reactions = 0

#..........................................................................#
# Mixture fraction setup
gas.TPY = T, 101000, oxidizer_inflow_composition

#elemental mass fraction of mt_element in oxidizer stream
Z_oxidizer = gas.elemental_mass_fraction(mixture_fraction_element)

gas.TPY = T, 101000, fuel_inflow_composition
#elemental mass fraction of mt_element in oxidizer stream
Z_fuel = gas.elemental_mass_fraction(mixture_fraction_element)





## ........................................................##
### TESTING ####


## Testing setup functions
def create_new_particle_list(length):
    particle_list = [None] * length
    for i,_ in enumerate(particle_list):
        particle_list[i] = pf.Particle(initial_composition, T, gas, Z_oxidizer,
                                       Z_fuel, mixture_fraction_element)
    return particle_list
# ................................................#

## Mixture fraction testing
p1 = pf.Particle(fuel_inflow_composition, 500, gas, Z_oxidizer,
                 Z_fuel, mixture_fraction_element)
p2 = pf.Particle(oxidizer_inflow_composition, T, gas, Z_oxidizer,
                Z_fuel, mixture_fraction_element)
pf.mix_particles( [(p1,p2)] )
# just baseline make sure it changes them
assert (p1.get_composition() != fuel_inflow_composition)
assert (p1.get_composition() != oxidizer_inflow_composition)


## Density testing
p1 = pf.Particle(fuel_inflow_composition, T, gas, Z_oxidizer,
                 Z_fuel, mixture_fraction_element)
p2 = pf.Particle(oxidizer_inflow_composition, 288.15, gas, Z_oxidizer,
                Z_fuel, mixture_fraction_element)
test_particle_list = [p1,p2]
print(pf.density(test_particle_list, gas))
# based on basic I.G.L. the densities should be related like this
assert(pf.density(test_particle_list, gas)[0] < pf.density(test_particle_list, gas)[1])
#now do it by hand and check the function
gas.TPY = p2.get_temperature(), 101000, p2.get_composition()
_,rho2 = gas.TD
assert(pf.density([p2],gas)[0] == rho2)


## Inflow/outlfow testing
replacement_particles = create_new_particle_list(50)
N_oxidizer_particles = round(0.5 * len(replacement_particles))
pf.inflow_outflow(replacement_particles, oxidizer_inflow_composition,
    fuel_inflow_composition, oxidizer_inflow_temperature,
    fuel_inflow_temperature, N_oxidizer_particles)
#make sure that all comps were changed to either fuel or oxidizer
for particle in replacement_particles:
    assert(
        particle.get_composition()['N2'] == oxidizer_inflow_composition['N2']
        or
        particle.get_composition()['N2'] == fuel_inflow_composition['N2']
        )


## Favre averged mixtre fraction testing
mixing_particle_list = create_new_particle_list(50)
original_mixture_fraction = pf.favre_averaged_mixture_fraction(
    mixing_particle_list, gas
)
mixing_pairs = [None] * round(len(mixing_particle_list) * 0.5)
for counter, mixing_pair in enumerate(mixing_pairs):
    mixing_pairs[counter] = (random.choice(mixing_particle_list),
        random.choice(mixing_particle_list))
pf.mix_particles(mixing_pairs)
new_mixture_fraction = pf.favre_averaged_mixture_fraction(
    mixing_particle_list, gas)
# famf should not change as all the particles are initially identical
assert (new_mixture_fraction == original_mixture_fraction)
#now try it with the inflow/outpflow and a mixing step
pf.inflow_outflow(mixing_particle_list, oxidizer_inflow_composition,
    fuel_inflow_composition, oxidizer_inflow_temperature,
    fuel_inflow_temperature, 0.5 * len(mixing_particle_list))
old_mixture_fraction = pf.favre_averaged_mixture_fraction(
    mixing_particle_list,gas)
mf1 = pf.mixture_fractions(mixing_particle_list,gas)
pf.mix_particles(mixing_pairs)
new_mixture_fraction = pf.favre_averaged_mixture_fraction(
    mixing_particle_list,gas)
mf2 = pf.mixture_fractions(mixing_particle_list,gas)

## Favre averaging test 2 a better test
p1 = pf.Particle(fuel_inflow_composition, T, gas, Z_oxidizer,
                 Z_fuel, mixture_fraction_element)
p2 =pf.Particle(oxidizer_inflow_composition, T, gas, Z_oxidizer,
                 Z_fuel, mixture_fraction_element)
p3 = pf.Particle(fuel_inflow_composition, T, gas, Z_oxidizer,
                 Z_fuel, mixture_fraction_element)
p4 = pf.Particle(oxidizer_inflow_composition, T, gas, Z_oxidizer,
                 Z_fuel, mixture_fraction_element)
#first do it by hand
gas.TPY = T, 101000, p1.get_composition()
_, rho1 = gas.TD
gas.TPY = T, 101000, p2.get_composition()
_, rho2 = gas.TD
mf1 = p1.get_mixture_fraction(gas)
mf2 = p2.get_mixture_fraction(gas)
famf1 = (0.5*(rho1*mf1 + rho2*mf2)) / (0.5* (rho1 + rho2))
#check the function
famf2 = pf.favre_averaged_mixture_fraction([p3,p4],gas)
#it would be great if these were the same
assert(famf1 == famf2)

density_list = pf.density(mixing_particle_list,gas)[20:30]
mixture_fraction_list = pf.mixture_fractions(mixing_particle_list, gas)[20:30]

#print (density_list, '\n', mixture_fraction_list)
famf = np.mean([rho*mf for rho,mf in zip(density_list, mixture_fraction_list)]) \
       / np.mean(density_list)

print(np.mean(mixture_fraction_list), famf)



### ..........................END OF TESTS......................##
print('\nTests Passed\n')
