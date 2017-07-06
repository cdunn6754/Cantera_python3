import numpy as np

### INPUTS

## Reactor Geometry
# Coal
r_i = 0.0013/2. #[m] radius of coal flow
A_coal = np.pi * r_i**2. #[m^2] area of coal flow

# Co-flow
A_coflow = 0.036**2.0 #[m^2] area of coflow
r_o = (A_coflow/np.pi + r_i**2)**0.5 # [m] outer radius of co-flow

## Flow Rates
Q_coalair = 5.0 #[SLPM]
Q_coflow = 10.0 #[SLPM]

## Densities
rho_ch4_standard = 0.664 #[kg/m^3] (IGL)
rho_air_standard = 1.20 #[kg/m^3]
rho_ch4 = 0.511 #[kg/m^3] #laramie pressure
rho_air = 0.924 #[kg/m^3]
coal_density = 1000.0 #[kg/m^3]

## Molecular wieghts
MW_ch4 = 16.04 #[g/mol]
MW_air = 28.97 #[g/mol]

## Coal flow mass fractions
Y_coal = 0.4
Y_coalair = 0.6

## Particle Sizes
particle_diameter = 70e-6 # [m] Diameter of coal particle

## OF wedge angle
wedge_angle = 10.0  # [deg]

## Equivelence ratio for co-flow (only works with methane)
phi = 0.8

## Equilibrium results from cantera for phi=0.8
T_ad = 1998.5 #[K]
rho_coflow_equil = 0.131 #[kg/m^3]



#................................Calculations...................................#
#...............................................................................#


### Coal Flow Calculations
#.........................

## Find the mass flow rates of the streams
m_dot_coalair = Q_coalair * (1.0/1000.) * (1./60.0) * rho_air_standard #[kg/s]

## Now find the velocity of the carrier gas
U_coalair = m_dot_coalair/(rho_air * A_coal) #[m/s]

## from mass fractions and mass flow rate of carrier gas find mdot of coal particles
m_dot_coal = m_dot_coalair * (0.4/0.6) #[kg/s]

## scale for axisymmetric 5 deg wedge
m_dot_coal_sim = m_dot_coal * (wedge_angle/360.0)


### Co-flow calculations
#.......................
 
# First come up with mass 
air_moles = 2. * 4.76 #[moles]
ch4_moles = phi #[moles]

# mole fration
X_ch4 = ch4_moles/(ch4_moles + air_moles)
X_air = air_moles/(ch4_moles + air_moles)
#mass fractions
Y_ch4 = (X_ch4 * MW_ch4)/(X_ch4 * MW_ch4 + X_air * MW_air)
Y_air = (X_air * MW_air)/(X_ch4 * MW_ch4 + X_air * MW_air)

## density of co-flow air at laramie pressure 
rho_coflow = rho_ch4 * Y_ch4 + rho_air * Y_air #[kg/m^3]
rho_coflow_standard = rho_ch4_standard * Y_ch4 + rho_air_standard * Y_air #[kg/m^3]

## Now find mass flow rate from Q_standard
m_dot_coflow = Q_coflow * (1.0/1000) * (1./60.0) * rho_coflow_standard #[kg/s]

print m_dot_coflow
## Find co-flow velocity (now using the equilibrium density) We model
## the coflow in OF as preburnt to it's equilibrium.
U_coflow = m_dot_coflow/(rho_coflow_equil * A_coflow) #[m/s]


### Random other stuff
#.....................

## Particle flow properties
particle_mass =  (4./3.) * np.pi * (particle_diameter/2.)**3.0 * coal_density #[kg]
particle_per_second = m_dot_coal_sim/particle_mass #[particle/second]


##............................OUTPUTS.....................................#
#...........................................................................#

print ('\ninflow velocity for co-flow: {} [m/s]'.format(U_coflow))
print ('inflow velocity for coal carrier air: {} [m/s]'.format(U_coalair))
print ('Coal particle m_dot (for 5 wedge): {} [kg/s]'.format(m_dot_coal_sim))
print ('Coal particle mass: {} [kg]'.format(particle_mass))
print ('Particles per second: {}\n'.format(particle_per_second))



