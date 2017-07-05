# Just wanted a function to read some text files (default to the local dirs I have here)
# but could be whatever. The funciton returns a string in the form that Cantera can
# use when your defining a composition (e.g. 'O2:0.3, N2:0.7' could be one for air).

# Pretty simple but I get sick of writing it in every single script. 

my_base_path = '/home/cdunn3/Python_2.7/Cantera/Volitile_composition/'
def read_compositions(file_name, base_dir = my_base_path):

    #open the files that describe the compositions
    fp = open(base_dir + file_name)

    # String form of the fuel composition to be returned
    composition = fp.read()

    # close the text file pointers
    fp.close()

    return composition

