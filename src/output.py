# Initialize an empty list to store all labels used in the plots
labs = []

# Define default labels for different parameters
elab = '$e$'
wlab = '$\omega$'
ilab = '$i$ (deg)'
klab = '$K$'
alab = '$a/R_\star$'

# Check for specific conditions and update labels accordingly
if is_ew:
    # If eccentricity and argument of periastron are combined
    elab = '$\sqrt{e} \sin \omega$'
    wlab = '$\sqrt{e} \cos \omega$'
if is_b_factor:
    # If impact parameter is used instead of inclination
    ilab = 'b'
if is_log_k:
    # If logarithm of the radial velocity semi-amplitude is used
    klab = '$\log_{10} K$'
if sample_stellar_density:
    # If stellar density is used instead of the semi-major axis to stellar radius ratio
    alab = '$\\rho_{\star}$'

# Generate labels for each planet's parameters
for i in range(nplanets):
    etiquetas = [
        f'$T0_{{{plabels[i]}}}$ [d]',  # Time of periastron passage
        f'$P_{{{plabels[i]}}}$ [d]',    # Orbital period
        f'{elab}$_{{{plabels[i]}}}$',   # Eccentricity or its variant
        f'{wlab}$_{{{plabels[i]}}}$',   # Argument of periastron or its variant
        f'{ilab}$_{{{plabels[i]}}}$',   # Inclination or impact parameter
        f'{alab}$_{{{plabels[i]}}}$' + '[${\\rm g\,cm^{-3}}$]',  # Semi-major axis ratio or stellar density
        f'{klab}$_{{{plabels[i]}}}$' + '[${\\rm km\,s^{-1}}$]'   # Radial velocity semi-amplitude or its log
    ]
    labs.append(etiquetas)

# Generate labels for each planet's radius
for i in range(nplanets):
    if nradius == 1:
        # Single radius label
        labs.append([f'$R_p/R_\star${plabels[i]}'])
    else:
        # Multiple radius labels for different bands
        for m in range(nradius):
            labs.append([f'$R_p/R_\star${plabels[i]}{bands[m]}'])

# Generate labels for limb darkening coefficients (LDCs)
for m in range(nbands):
    labs.append([f'$q_1${bands[m]}', f'$q_2${bands[m]}'])

# Add labels for radial velocity (RV) instruments
labs.append(telescopes_labels)

# Add labels for RV jitter and transit jitter
for i in range(n_jrv):
    labs.append([f'$\\sigma_j$, {rv_jitter_labels[i]}'])
for i in range(n_jtr):
    labs.append([f'$\\sigma_j$, {bands[i]}'])

# Add labels for linear and quadratic trends
labs.append(['Linear trend'])
labs.append(['Quadratic trend'])

# Add labels for Gaussian Process RV and transit model parameters
labs.append(gprv_labels)
labs.append(ktr_labels)

# Concatenate all the generated labels into a single array
labels = np.concatenate(labs)


# ===========================================================
#              plot chains
# ===========================================================

#Load plot routines
exec(open('src/plotting_functions.py').read())


if (len(plot_parameters) < 2):

    # plot_parameters stores the flag of each fitted parameter
    plot_parameters = []

    #We only plot sampled parameters
    for o in range(len(prior_flags)):
        if (prior_flags[o] != 'f'):
            plot_parameters.append(o)

#Let us compute a numpy matrix to create the plots, this works better with the functions
df_np = np.array(df[sampled_parameters]).T

if is_plot_chains:
    create_chains_plot(df.chain_number,np.array(df.log_likelihood),df_np,labels,plot_parameters)

if is_plot_posterior:
    create_plot_posterior(df_np, labels,cbars='red', nb=50, num=plot_parameters)

if is_plot_correlations:
    create_plot_correlation(df_np, labels, col='blue', num=plot_parameters)

#if is_corner_plot:
#    create_corner_plot(df_np, labels,maxloglike,plot_parameters)

#Compute the median of all the parameters to be used to compute the plots
df_medians = df.median()

#To plot random samples of the posterior
n_samples = 100
df_samples = df.sample(n=n_samples)

 # Transit plots
if total_tr_fit:
    exec(open('src/plot_tr.py').read())

# RV plots
if total_rv_fit:
    exec(open('src/plot_rv.py').read())

sys.exit('Exit at the plot stage succesfully')

