#!/usr/bin/env python

# Dakota will execute this script as
#  The command line arguments will be extracted by dakota.interfacing automatically.

# necessary python modules
import dakota.interfacing as di

# ----------------------------
# Parse Dakota parameters file
# ----------------------------

params, results = di.read_parameters_file()

# -------------------------------
# Convert and send to application
# -------------------------------

# set up the data structures 
# for this simple example, put all the variables into a single hardwired array
# The asv has to be mapped back into an integer
continuous_vars = [ params['delta_r2'],params['delta_r3'],params['r_diameter'],params['r_height']]
active_set_vector = 0
if results["obj_fn"].asv.function:
    active_set_vector += 1


# Alternatively, the ASV can be accessed by index in
# function, gradient, hessian order
#for i, bit in enumerate(results["obj_fn"].asv):
#    if bit:
#        active_set_vector += 1 << i

# set a dictionary for passing to via Python kwargs
solstice_params = {}
solstice_params['cv'] = continuous_vars
solstice_params['asv'] = [active_set_vector]
solstice_params['functions'] = 1

# execute the analysis as a separate Python module
print "Running dakota..."
from one_key_co_optimisation import one_key_dakota

solstice_results = one_key_dakota(**solstice_params)
print "complete."

# ----------------------------
# Return the results to Dakota
# ----------------------------

# Insert functions from rosen into results
# Results.responses() yields the Response objects.
for i, r in enumerate(results.responses()):
    if r.asv.function:
        r.function = solstice_results['fns'][i]

results.write()
