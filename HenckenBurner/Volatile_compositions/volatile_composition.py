import numpy as np


# some molar weights
M_CO = 28.01
M_CO2 = 44.01
M_CH4 = 16.04
M_H2 = 2.01
M_C2H4 = 28.05
M_C6H6 = 78.11 # tar
M_H2S = 34.08
M_N2 = 28.02
M_C = 12.01
M_H = 1.01
M_O = 16.0
M_N = 14.01
M_S = 32.1


## P and U analyses are from Belmont PRB email
# Proximate analysis daf # mass basis
p_volatile_daf = 0.4834
p_fixedC_daf = 1.0 - p_volatile_daf
# total volatile yield = tvy*promimate yield
Q = 1.2 # guess
p_volatile = p_volatile_daf*Q
p_fixedC = 1.0 - p_volatile


# ultimate analysis daf # mass basis
u_O = 0.1389
u_H = 0.0423
u_N = 0.0124
u_S = 0.0000 ## there actaully is u_S = 0.0050 but mech can't handle it
u_C = 1-(u_O+u_H+u_N+u_S)


# calculate carbon in the volatiles: C from ultimate - fixed carbon C
v_C = u_C - p_fixedC
# calculate mass fractions in volatiles
sum = u_O + u_H + u_N + u_S + v_C
v_O = u_O*1/sum
v_H = u_H*1/sum
v_N = u_N*1/sum
v_S = u_S*1/sum
v_C = v_C*1/sum
# calculate elemental mol fractions in volatiles
sum = v_C/M_C + v_H/M_H + v_O/M_O + v_N/M_N + v_S/M_S
nC = v_C/M_C / sum
nH = v_H/M_H / sum
nO = v_O/M_O / sum
nN = v_N/M_N / sum
nS = v_S/M_S / sum

elemental_mol_fractions = [nH, nO, nC, nN] 

# O_splitter = 0.01
# C_splitter = 0.5

# # molar fractions of volatiles
# X_N2   = 0.5*nN
# X_CO  = O_splitter*nO
# X_CO2  = 0.5*(1. - O_splitter)*nO
# X_C2H4  = 0.0 #0.5*(C_splitter*nC - X_CO2 - X_CO)
# X_CH4 = C_splitter*nC - 2*X_C2H4 - X_CO2 - X_CO
# X_H2   = nH - 4.*X_C2H4 - 4.*X_CH4
# X_C6H6 = 
# sum = X_N2 + X_CO2 + X_C2H4 + X_H2 + X_CO + X_CH4
# X_N2   = X_N2/sum
# X_CO2  = X_CO2/sum
# X_C2H4  = X_C2H4/sum
# X_CH4 = X_CH4/sum
# X_H2   = X_H2/sum
# X_CO = X_CO/sum

# Good mixture for a lot of benzene
# this fuel has LHV of 27.95 [Mj/kg]
# xiCO = 0.0
# xiTar = 0.1
# xiC2H4 = 0.0899
# fuel_name = 'heavy'

#No Tar
#LIGHT
xiCO = 0.95
xiTar = 0.07
xiC2H4 = 0.099
fuel_name = 'light'

# molar fractions of volatiles
X_N2   = 0.5*nN
X_CO   = xiCO*nO
X_CO2  = 0.5*(1-xiCO)*nO
X_CH4  = (1 - 2*xiC2H4 - 6*xiTar)*nC - 0.5*(xiCO + 1)*nO
X_H2   = 0.5*nH - 2*(1 - xiC2H4 - 4.5*xiTar)*nC + (xiCO+1)*nO - nS
X_C2H4 = xiC2H4*nC
X_C6H6  = xiTar*nC
sum = X_N2 + X_CO + X_CO2 + X_CH4 + X_H2 + X_C2H4 + X_C6H6
X_N2   = X_N2/sum
X_CO   = X_CO/sum
X_CO2  = X_CO2/sum
X_CH4  = X_CH4/sum
X_H2   = X_H2/sum
X_C2H4 = X_C2H4/sum
X_C6H6 = X_C6H6/sum

mol_fractions = {'N2':X_N2, 'CO2': X_CO2, 'CH4':X_CH4,
                 'H2':X_H2, 'CO':X_CO, 'C2H4':X_C2H4,
                 'C6H6': X_C6H6}

print("The mol fractions: ")
for specie in mol_fractions:
    print( specie, mol_fractions[specie])

# here verify that the elemental mol fractions of the volatiles
# are consistent with inputs
elemental = np.zeros(4)
elemental[0] = 4.0*X_C2H4 + 2.0*X_H2 + 4.*X_CH4 + 6*X_C6H6#H
elemental[1] = 2.0*X_CO2 + X_CO #O
elemental[2] = 1.0*X_CO2 + 2.0*X_C2H4 + X_CO  + 6*X_C6H6 + X_CH4#C
elemental[3] = 2.0*X_N2 #N

print ("\nDifference: ", elemental_mol_fractions  - elemental/np.sum(elemental))


# mass fractions of volatiles
sum = X_N2*M_N2 + X_CO*M_CO + X_CO2*M_CO2 + X_H2*M_H2 + X_C2H4*M_C2H4 + X_CH4*M_CH4 + X_C6H6*M_C6H6
Y_N2   = X_N2*M_N2/sum
Y_CO   = X_CO*M_CO/sum
Y_CO2  = X_CO2*M_CO2/sum
Y_H2   = X_H2*M_H2/sum
Y_C2H4 = X_C2H4*M_C2H4/sum
Y_CH4 = X_CH4*M_CH4/sum
Y_C6H6 = X_C6H6*M_C6H6/sum

species_mass_fractions = {'N2':Y_N2, 'CO2': Y_CO2, 'CH4':Y_CH4,
                          'H2':Y_H2, 'CO': Y_CO, 'C2H4': Y_C2H4,
                          'C6H6': Y_C6H6}

print("\nThe mass fractions: ")
for specie in species_mass_fractions:
    print( specie, species_mass_fractions[specie])

# Calculate the lower heating value of the gas mixture: should be equal to
# the "measured" LHV of the volatiles
LHV_coal = 29.28    # [Mj/kg] heating value of coal -> from Belmont analysis email
LHV_char = 32.9    # heating value of char
LHV_measured_volatiles = (LHV_coal - LHV_char*p_fixedC)/p_volatile

# We need the LHV's of the different species
LHV_H2S = 0
LHV_N2 = 0
LHV_CO = 10.112
LHV_CO2 = 0
LHV_CH4 = 50.009
LHV_H2 = 119.96
LHV_C2H4 = 47.195
LHV_C6H6 = 40.5

LHV_vols = Y_N2*LHV_N2 + Y_CO*LHV_CO + Y_CO2*LHV_CO2 +  Y_CH4*LHV_CH4 +  \
           Y_H2*LHV_H2 + Y_C2H4*LHV_C2H4 + Y_C6H6*LHV_C6H6

print ('\nMeasured Volatile LHV: {}'.format(LHV_measured_volatiles))
print ('LHV for this Mixture: {}'.format(LHV_vols))
#print('This should be positive or 0.0:  {}'.format(np.min(mol_fractions)))
exit()

# # write out a string of molar fractions that can be directly imported to cantera scripts
# fp = open('%s_fuel.txt' %fuel_name,'w+')
# fp.write('N2:%f,CO:%f,CO2:%f,CH4:%f,H2:%f,C2H4:%f,C6H6:%f'  \
#         % ( X_N2, X_CO, X_CO2, X_CH4, X_H2, X_C2H4, X_C6H6))
# fp.close()
