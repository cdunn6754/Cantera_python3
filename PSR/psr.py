"""
A PSR combustor. Two separate streams - one coal volatiles and the other air, both at
300 K and 1 atm flow into an adiabatic combustor where they mix and burn.

We are interested in the steady-state burning solution. Since at 300 K no
reaction will occur between fuel and air, we need to use an 'igniter' to
initiate the chemistry. A simple igniter is a pulsed flow of atomic hydrogen.
After the igniter is turned off, the system approaches the steady burning
solution.
"""

# let me get at my utility functions
import sys 
import os
sys.path.append(os.path.abspath("/home/cdunn3/Python_3.6.1/Cantera/utilities"))

import math
import csv
import matplotlib.pyplot as plt
import cantera as ct
import numpy as np
import CompositionReader as cr

#..........................................................................#
## setup stuff

# user inputs
mechanism = 'reduced'
fuel_name = 'heavy'

# Get the specified fuel compositions
fuel_dict = {'light':'light_fuel.txt',
             'heavy':'heavy_fuel.txt'
         }
composition = cr.read_compositions(fuel_dict[fuel_name])

# define the mechanism and gas
mech_dir = 'mechanisms/'
mechanism_dict = {'reduced':'r_creck_52.cti',
                  'full':'POLIMI_TOT_1412.cti'
              }
gas = ct.Solution(mech_dir +  mechanism_dict[mechanism])

# ....................................................................#
## run the simulation

# create a reservoir for the fuel inlet, and set to pure methane.
gas.TPX = 300.0, ct.one_atm, composition
fuel_in = ct.Reservoir(gas)
fuel_mw = gas.mean_molecular_weight

# use predefined function Air() for the air inlet
air = ct.Solution('air.xml')
air_in = ct.Reservoir(air)
air_mw = air.mean_molecular_weight

# to ignite the fuel/air mixture, we'll introduce a pulse of radicals. The
# steady-state behavior is independent of how we do this, so we'll just use a
# stream of pure atomic hydrogen.
gas.TPX = 300.0, ct.one_atm, 'H:1.0'
igniter = ct.Reservoir(gas)

# create the combustor, and fill it in initially with N2
gas.TPX = 300.0, ct.one_atm, 'N2:1.0'
combustor = ct.IdealGasReactor(gas)
combustor.volume = 1.0

# create a reservoir for the exhaust
exhaust = ct.Reservoir(gas)

# adjust the ratio coming in
Z  = 0.15

# compute fuel and air mass flow rates (light crashes at 2.5)
total_mdot = lambda t:  1 + (2.*t)**2.0
air_mdot = lambda t: total_mdot(t) * (1.0 -Z)
fuel_mdot = lambda t: total_mdot(t) * Z

# create and install the mass flow controllers. Controllers m1 and m2 provide
# constant mass flow rates, and m3 provides a short Gaussian pulse only to
# ignite the mixture
m1 = ct.MassFlowController(fuel_in, combustor, mdot=fuel_mdot)

# note that this connects two reactors with different reaction mechanisms and
# different numbers of species. Downstream and upstream species are matched by
# name.
m2 = ct.MassFlowController(air_in, combustor, mdot=air_mdot)

# The igniter will use a Gaussian time-dependent mass flow rate.
fwhm = 0.2
amplitude = 0.01
t0 = 0.5
igniter_mdot = lambda t: amplitude * math.exp(-(t-t0)**2 * 4 * math.log(2) / fwhm**2)
m3 = ct.MassFlowController(igniter, combustor, mdot=igniter_mdot)

# put a valve on the exhaust line to regulate the pressure
v = ct.Valve(combustor, exhaust, K=1.0)

# the simulation only contains one reactor
sim = ct.ReactorNet([combustor])

# take single steps to 6 s, writing the results to a CSV file for later
# plotting.
tfinal = 20.0
dt = 0.0001
times = np.arange(dt,tfinal,dt)

# data output arrays
states = ct.SolutionArray(gas, extra=['t', 'tres', 'fuel_mdot', 'Z'])
for i, time in enumerate(times):
    sim.advance(time)
    tres = combustor.mass/(m1.mdot(time) + m2.mdot(time))    #v.mdot(time)
    states.append(gas.state, t=time, tres=tres, fuel_mdot=m1.mdot(time), Z=Z)

csv_dir = 'csv_files/'
states.write_csv(csv_dir + fuel_name + '_' + mechanism + '_states.csv')
