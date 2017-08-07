import pickle
import particle_functions as pf

particle_list = pickle.load(open('saved_particles.p', 'rb'))

print(pf.average_properties(particle_list))
