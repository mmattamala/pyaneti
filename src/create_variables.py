import pickle

# Load the dictionary from the pickle file
fpickle = f'{outdir}/posterior_parameters.pkl'
with open(fpickle, 'rb') as pickle_file:
    loaded_data = pickle.load(pickle_file)

# Assign dictionary elements to variables
for key, value in loaded_data.items():
    locals()[key] = value
