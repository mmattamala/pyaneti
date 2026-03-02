#Let us define a likelihood that can be used in python calling the fortran routines
def pyaneti_likelihood(pars,rv_time,rv_vals,rv_errs,lc_time,lc_flux,lc_errs,
        rvlab,jrvlab,dim_lab,trlab,jtrlab,tff,flags,kernels,model_int,model_double):

        nll = pti.get_loglike(rv_time,rv_vals,lc_time,lc_flux,rv_errs,lc_errs,
        rvlab,jrvlab,dim_lab,trlab,jtrlab,tff,flags,kernels,pars,model_int,model_double)[0]

        if not np.isfinite(nll):
            return -np.finfo(float).max

        return nll

#Let us define a likelihood that can be used in python calling the fortran routines
def pyaneti_likelihood_parallel(pars_to_sample):

        #Prepare the data each time the function runs, this helps with parallel running
        #exec(open('src/create_variables.py').read())

        #Not all parameters that pyaneti needs are sampled, so let us transform the pars_to_sample
        #parameter into pars that has to be passed to pyaneti
        pars = np.zeros(shape=(len(prior_flags)))
        values = prior_vals[::2]
        ifixed = prior_flags == 'f'
        pars[ifixed] = values[ifixed]
        pars[~ifixed] = pars_to_sample

        nll = pti.get_loglike(rv_time,rv_vals,lc_time,lc_flux,rv_errs,lc_errs,
        tlab,jrvlab,dim_lab,trlab,jtrlab,total_fit_flag,flags,kernels,pars,model_int,model_double)[0]

        return nll

def pyaneti_priors(pars,prior_flags,prior_vals):
        return pti.get_priors(prior_flags,prior_vals,pars)

def pyaneti_posterior(pars_to_sample,rv_time,rv_vals,rv_errs,lc_time,lc_flux,lc_errs,
        prior_flags, prior_vals,
        rvlab,jrvlab,dim_lab,trlab,jtrlab,tff,flags,kernels,model_int,model_double):

        #Not all parameters that pyaneti needs are sampled, so let us transform the pars_to_sample
        #parameter into pars that has to be passed to pyaneti
        pars = np.zeros(shape=(len(prior_flags)))
        values = np.array(prior_vals[::2])
        ifixed = prior_flags == 'f'
        pars[ifixed] = values[ifixed]
        pars[~ifixed] = np.array(pars_to_sample)

        nll = pyaneti_likelihood(pars,rv_time,rv_vals,rv_errs,lc_time,lc_flux,lc_errs,
        rvlab,jrvlab,dim_lab,trlab,jtrlab,tff,flags,kernels,model_int,model_double)

        if nll != nll:
            return -np.finfo(float).max

        priors = pti.get_priors(prior_flags,prior_vals,pars)
        if any(priors < 1e-6):
            return -np.finfo(float).max

        log_priors = np.sum(np.log(priors))


        return log_priors + nll

def pyaneti_posterior_parallel(pars_to_sample):

        #Prepare the data each time the function runs, this helps with parallel running
        #exec(open('src/create_variables.py').read())

        posterior = pyaneti_posterior(pars_to_sample,rv_time,rv_vals,rv_errs,lc_time,lc_flux,lc_errs,
        prior_flags, prior_vals,tlab,jrvlab,dim_lab,trlab,jtrlab,total_fit_flag,flags,kernels,model_int,model_double)

        return posterior


#This function transform priors for the nested samplig routienes that sample the unit multidimensional cube
def transform_priors(cube):
        import scipy.stats

        #Prepare the data each time the function runs, this helps with parallel running
        #exec(open('src/create_variables.py').read())

        #Let us do the transformation for the parameters that we are sampling only
        isampled = prior_flags != 'f'
        sampled_flags = prior_flags[isampled]
        tpriors = np.zeros(shape=(len(sampled_flags)))
        j = 0
        for i,flag in enumerate(prior_flags):
            if flag == 'u':
                tpriors[j] = cube[j] * (prior_vals[2*i+1] - prior_vals[2*i]) + prior_vals[2*i]
                j += 1
            elif flag == 'g':
                gaussian = scipy.stats.norm(prior_vals[2*i],prior_vals[2*i+1])
                tpriors[j] = gaussian.ppf(cube[j])
                j += 1
            elif flag == 'j': #Jeffrey's prior
                lower_bound = prior_vals[2 * i]
                upper_bound = prior_vals[2 * i + 1]
                tpriors[j] = lower_bound * (upper_bound / lower_bound)**cube[j]
                j += 1
            elif flag == 'm':  # Modified Jeffrey's prior
                lx = prior_vals[2 * i]
                rx = prior_vals[2 * i + 1]
                # Transform the cube value using the inverse CDF method
                log_term = np.log((lx + rx) / lx)
                tpriors[j] = lx * ((1.0 - np.exp(-cube[j] * log_term)) / np.exp(-cube[j] * log_term))
                j += 1
            elif flag == 'b':  # Beta distribution
                alpha = prior_vals[2 * i]
                beta = prior_vals[2 * i + 1]
                tpriors[j] = scipy.stats.beta.ppf(cube[j], alpha, beta)
                j += 1
            elif flag == 'f':
                continue
            else:
                raise ValueError(f"Nested sampling prior  tranformation for '{flag}' flag not defined yet")

        return tpriors


if (method == 'mcmc'):

        # Ensure nwalkers is divisible by 2
        if (nwalkers % 2 != 0):
            nwalkers = nwalkers + 1

        pti.mcmc_stretch_move(
            rv_time, rv_vals, lc_time, lc_flux, rv_errs, lc_errs,
            tlab, jrvlab, dim_lab, trlab, jtrlab,
            flags, total_fit_flag, prior_flags, prior_vals,
            kernels, model_int,
            model_double,
            nwalkers, maxi, thin_factor, nconv)

        #Change fortran file to a lighter csv file
        df = pd.read_csv('all_data.dat',header=None,sep='\s+')
        df.to_csv('all_data.dat', header=False,index=False,sep=',')

elif method in ['zeus','emcee']:


        if method == 'zeus': import zeus as the_sampler
        if method == 'emcee': import emcee as the_sampler
        import multiprocessing
        from multiprocessing import Pool
        import os
        os.system("export OMP_NUM_THREADS=1")

        #Let us extract the indices of the parameters that we are sampl
        sampled_indices = prior_flags != 'f'
        #Dimension of the sampled space
        ndim = len(prior_flags[sampled_indices])

        print(f"Running MCMC sampling using {method}")
        print(f"Doing {2*thin_factor*nconv} steps, sampling with {nwalkers} chains")
        print(f"Sampling {ndim} parameters")
        if not is_parallel_run: print("Sequential run")
        if is_parallel_run: print("Parallel run")

        #Number os steps
        nsteps = 2*thin_factor*nconv
        #Initiate the chains using the create_chains function of pyaneti
        #Let us do a while cycle to ensure that all intial chains have valid solutions for the posterior
        pars_to_sample = np.zeros(shape=(nwalkers,ndim))
        i = 0
        while i < nwalkers:
            k = 0
            for j in range(len(prior_flags)):
                if prior_flags[j] != 'f':
                    pars_to_sample[i,k] = pti.create_chains(prior_flags[j],prior_vals[2*j:2*j+2])[0]
                    k += 1
            #Estimate the posterior we get from this chain, if it is ridiculously small, we go for otherev
            lp = pyaneti_posterior(pars_to_sample[i,:],rv_time,rv_vals,rv_errs,lc_time,lc_flux,lc_errs,
            prior_flags, prior_vals,
            tlab,jrvlab,dim_lab,trlab,jtrlab,total_fit_flag,flags,kernels,model_int,model_double)
            if lp > -1e100:
                i += 1

        if is_parallel_run:
            with Pool() as pool:
                    sampler = the_sampler.EnsembleSampler(nwalkers, ndim, pyaneti_posterior_parallel,
                    pool=pool)
                    sampler.run_mcmc(pars_to_sample, nsteps, progress=True)
        else:
            sampler = the_sampler.EnsembleSampler(nwalkers, ndim, pyaneti_posterior,args=
            (rv_time,rv_vals,rv_errs,lc_time,lc_flux,lc_errs,
            prior_flags, prior_vals,
            tlab,jrvlab,dim_lab,trlab,jtrlab,total_fit_flag,flags,kernels,model_int,model_double)
            )
            sampler.run_mcmc(pars_to_sample, nsteps, progress=True)

        log_prob = sampler.get_log_prob()

        if method == "zeus": chains = sampler.chain
        if method == "emcee": chains = sampler.get_chain()

        #Let's create the file that will be used for pyaneti outputs
        table = []
        for i in range(thin_factor*nconv,nsteps,thin_factor):
            for j in range(nwalkers):
                sampled_pars = np.zeros(shape=(len(prior_flags)))
                l = 0
                for k,flag in enumerate(prior_flags):
                    if flag == 'f': #the parameter is fixed, let us assing the value
                        sampled_pars[k] = prior_vals[2*k]
                    else: #the parameter is being sampled, let us give it the value that was sampled
                        sampled_pars[k] =chains[i,j,l]
                        l += 1
                table.append([i,j,log_prob[i,j],*sampled_pars])

        df = pd.DataFrame(table)
        df.to_csv('all_data.dat', header=False,index=False,sep=',')


elif method == 'ultranest':

        import ultranest

        sampled_indices = prior_flags != 'f'
        #Extract names of the sampled parameters
        param_names = list(labs[sampled_indices])
        sampler = ultranest.ReactiveNestedSampler(param_names,pyaneti_likelihood_parallel, transform_priors)
        result = sampler.run()
        sampler.print_results()

        sys.exit(f"{method} method ended succesfully. {method} run does not generate plots.")

elif (method == 'dynesty'):

        import dynesty
        from multiprocessing import Pool

        #Let us extract the indices of the parameters that we are sampl
        sampled_indices = prior_flags != 'f'
        #Dimension of the sampled space
        ndim = len(prior_flags[sampled_indices])

        dsampler = dynesty.NestedSampler(pyaneti_likelihood_parallel, transform_priors, ndim)
        dsampler.run_nested(print_progress=True)


        # Get the results
        results = dsampler.results

        # Print a summary of the results
        results.summary()

        sys.exit(f"{method} method ended succesfully. {method} run does not generate plots.")


elif (method == 'plot'):
        print('I will only print the values and generate the plots')
        newfile = outdir+'/'+star+'_all_data.dat'
        if not os.path.exists(newfile):
            print(f"The file {newfile} for {star} does not exist.")
            print("Create it using one of the valid MCMC sampling methods")
            print('method = mcmc      -> Run the MCMC code using the stretch move algorithm.')
            print('method = zeus      -> Run the MCMC sampling using the Zeus sampler.')
            print('method = emcee     -> Run the MCMC sampling using the emcee sampler.')
            sys.exit("ERROR!")

else:
    print('You did not choose a valid method!')
    print('Please choose from the following methods:')
    print('method = mcmc      -> Run the MCMC code using the stretch move algorithm.')
    print('method = zeus      -> Run the MCMC sampling using the Zeus sampler.')
    print('method = emcee     -> Run the MCMC sampling using the emcee sampler.')
    print('method = ultranest -> Run the Nested Sampling using the UltraNest sampler.')
    print('method = dynesty   -> Run the Nested Sampling using the Dynesty sampler.')
    print('method = plot      -> Generate plots from a previous run.')
    sys.exit('Please choose your favorite method and rerun the script.')



#Prepare the data file and labels for the upcoming part of the excecution
#Printing the values, compute statistical ranges, create plots, etc.
newfile = outdir+'/'+star+'_all_data.dat'
if (os.path.isfile('all_data.dat')):
        os.rename('all_data.dat', newfile)

        labels = np.concatenate([['i'],['chain_number'], ['log_likelihood'], labs])
        labels = [item + ',' for item in labels]
        labels = [item.strip() for item in labels]
        line_prepender(newfile, labels)
