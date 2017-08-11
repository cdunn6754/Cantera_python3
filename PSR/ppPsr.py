# Clint Dunn, cdunn6754@gmail.com, 07-5-17
# written in Python 3.6.1 on cantera 2.3.0


##  Description
# Script for plotting based on the states.csv files that psr.py produces.
# I would like to seperate plotting from the actual simulations which can be 
# quite time consuming. So I just run the simulations once and then 
# can mess around in this script with plotting stuff

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

## read the csv files to pandas dataframes
csv_dir = 'csv_files/'
light_reduced_states= pd.read_csv(csv_dir + "light_drm_states.csv")
light_full_states = pd.read_csv(csv_dir + "light_full_states.csv")
heavy_reduced_states= pd.read_csv(csv_dir + "heavy_drm_states.csv")
heavy_full_states= pd.read_csv(csv_dir + "heavy_full_states.csv")

## get data numpy arrays
#time
times = light_full_states.loc[:,'t'].as_matrix() #the times will be the same for all
long_times = heavy_reduced_states.loc[:,'t'].as_matrix() # times for the longer
# drm case

#temperatures
light_full_temps = light_full_states.loc[:,'T'].as_matrix()
light_reduced_temps = light_reduced_states.loc[:,'T'].as_matrix()
heavy_full_temps = heavy_full_states.loc[:,'T'].as_matrix()
heavy_reduced_temps = heavy_reduced_states.loc[:,'T'].as_matrix()

#residence times
heavy_reduced_tres = heavy_reduced_states.loc[:,'tres'].as_matrix()

## plotting
# all combined in a temp vs. time plot
plt.figure(0)
plt.plot(times,light_reduced_temps, label='light fuel; reduced mechanism')
plt.plot(times,light_full_temps, label='light fuel; full mechanism')
plt.plot(long_times,heavy_reduced_temps, label='heavy fuel; reduced mechanism')
plt.plot(times,heavy_full_temps, label='heavy fuel; full mechanism')
plt.title('Extinction behavior comparison of mechanism and fuel composition')
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.legend()

# residence times
plt.figure(1)
plt.plot(long_times,heavy_reduced_tres, label='heavy fuel; reduced mechanism')
# plt.plot(times,light_full_temps, label='light fuel; full mechanism')
# plt.plot(times,heavy_reduced_temps, label='heavy fuel; reduced mechanism')
# plt.plot(times,heavy_full_temps, label='heavy fuel; full mechanism')
plt.title('Extinction behavior comparison of mechanism and fuel composition')
plt.xlabel('Time [s]')
plt.ylabel('Residence Times [s]')
plt.ylim([0,0.05])
plt.legend()



plt.show()

