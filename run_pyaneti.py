# -----------------------------------------------------------
#                         pyaneti.py
#                     Main pyaneti file
#                   Barragan O, March 2016
# -----------------------------------------------------------

# Load libraries
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import pyaneti as pti  # FORTRAN module

# -------------------------------------------------------------
#                   INITIALIZATION
# -------------------------------------------------------------

# Load the input file
# You have to run the program as ./pyaneti star_name
star = str(sys.argv[1])

# Create path to the input_fit.py file
inf_name = f'inpy/{star}/input_fit.py'

# Did you create an input_fit.py file?
if (not os.path.isfile(inf_name)):
    print('You have not created', inf_name)
    sys.exit()

# Read the file with all the python functions
exec(open('src/todo-py.py').read())

# Read the file with the default values
exec(open('src/default.py').read())

# Read input file
exec(open(inf_name).read())

# Create ouput directory
outdir = f'{outdir}{star}_out'
if not os.path.exists(outdir):
    os.makedirs(outdir)

#Read the data and prepare all the variables
exec(open('src/prepare_data.py').read())

# -------------------------------------------------------------
#                   SAMPLING ROUTINES
# -------------------------------------------------------------

exec(open('src/samplers.py').read())

# -------------------------------------------------------------
#               PRINT AND PLOT ROUTINES
# -------------------------------------------------------------

# Print the values
exec(open('src/print_values.py').read())

exec(open('src/output.py').read())

# -------------------------------------------------------------
#             	 END pyaneti.py FILE
# -------------------------------------------------------------
