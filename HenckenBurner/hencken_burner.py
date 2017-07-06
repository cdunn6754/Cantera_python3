"""
Adiabatic flame temperature and equilibrium composition for a fuel/air mixture
as a function of equivalence ratio, including formation of solid carbon.
"""

import cantera as ct
import numpy as np
import sys
import csv
import matplotlib.pyplot as plt

##############################################################################
# Edit these parameters to change the initial temperature, the pressure, and
# the phases in the mixture.

T = 295 #room temp
P = 78000 #laramie atm

species = {S.name: S for S in ct.Species.listFromFile('r_creck_52.cti')}

# phases
# complete_species = [species[S] for S in ('CH4','O2','N2','CO2','H2O')]
# gas = ct.Solution(thermo='IdealGas', species=complete_species)
gas = ct.Solution('r_creck_52.cti')
carbon = ct.Solution('graphite.xml')

# the phases that will be included in the calculation, and their initial moles
mix_phases = [(gas, 1.0), (carbon, 0.0)]

# gaseous fuel species
fuel_species = 'CH4'


##############################################################################

mix = ct.Mixture(mix_phases)

# create some arrays to hold the data
phi = 0.8

# find fuel, nitrogen, and oxygen indices
ifuel = gas.species_index(fuel_species)
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

stoich_O2 = gas.n_atoms(fuel_species,'C') + 0.25*gas.n_atoms(fuel_species,'H')

X = np.zeros(gas.n_species)
X[ifuel] = phi
X[io2] = stoich_O2
X[in2] = stoich_O2*3.76

# set the gas state
gas.TPX = T, P, X

# create a mixture of 1 mole of gas, and 0 moles of solid carbon.
mix = ct.Mixture(mix_phases)
mix.T = T
mix.P = P
print ('\nIncoming gas mole fractions:')
print (gas.X)
print ('\n')

# equilibrate the mixture adiabatically at constant P
mix.equilibrate('HP', solver='gibbs', max_steps=1000)

t_ad = mix.T
print ('The adiabatic flame temperature: {:0.1f}'.format(t_ad))
print ('The equilibrium density: {:0.3f}'.format(gas.density))


# sort through the species mass fractions to find those which
# exceed the threshold eps, then store them along with their names in
# dictionary "relevant_species", I think eps==1e-6 is good
eps = 1e-4
relevant_species = {}
mass_fraction_sum = 0.
for i in range(len(gas.species_names)):
    if gas.Y[i] >= eps:
        current_name = gas.species_names[i]
        current_mass_fraction = gas.Y[i]
        mass_fraction_sum = mass_fraction_sum + current_mass_fraction
        relevant_species[current_name] = current_mass_fraction

# renormalize the remaining mass fractions
for key in relevant_species:
    relevant_species[key] = relevant_species[key]/mass_fraction_sum

# print them out
print ('\nNumber of relevant species: {}'.format(len(relevant_species)))
for key in relevant_species:
    print (key,':\t' ,relevant_species[key])
        


    
