# Let us do the plots here
from matplotlib import gridspec
from matplotlib.colors import LogNorm

if (is_seaborn_plot):
    import seaborn as sns
    sns.set(style=seaborn_style)
    sns.set_color_codes(seaborn_palette)
    sns.set_context(seaborn_context)

fsx = figure_size_x
fsy = figure_size_y
fos = font_size_label


def create_chains_plot(chain_number, posterior, params, labels, plot_parameters):

    # Reduced figure size to 1/4th of the original size
    fname = f'{outdir}/{star}_chains.pdf'
    print('Creating ', fname)

    #Extract the chains available in chain_number
    chains = list(set(chain_number))

    # Original figsize divided by 2 for width and height (results in 1/4 of area)
    plt.figure(1, figsize=(fsx, len(plot_parameters) * fsy / 4))

    # Setup grid for subplots
    gs = gridspec.GridSpec(nrows=len(plot_parameters) + 1, ncols=1)
    gs.update(hspace=0.0)  # Small hspace for compactness

    # Set font sizes for labels and tick numbers
    label_font_size = 8  # Font size for labels
    tick_font_size = 6   # Font size for tick numbers

    # Plot the posterior
    ax0 = plt.subplot(gs[0], rasterized=is_rasterized)
    ax0.set_ylabel('log_likelihood', fontsize=label_font_size)
    n_iters = sum(0 == chain_number)
    for i in chains:
        indices = chain_number == i
        ax0.plot(posterior[indices], alpha=0.3, lw=0.3)
    ax0.set_xticks([])  # No x-ticks on the top plot
    ax0.tick_params(axis='y', which='both', labelsize=tick_font_size)
    ax0.set_xlim(0, n_iters-1)
    plt.ticklabel_format(useOffset=False, axis='both')

    # Plot the parameters
    for n, param in enumerate(plot_parameters, start=1):
        ax = plt.subplot(gs[n], rasterized=is_rasterized)
        ax.set_ylabel(labels[param], fontsize=label_font_size)  # Smaller font size for ylabel
        for i in chains:
            indices = chain_number == i
            ax.plot(params[param][indices], alpha=0.3, lw=0.3)
        if param != plot_parameters[-1]:
            ax.set_xticks([])  # Skip x-ticks for intermediate plots
        else:
            ax.set_xlabel('iteration', fontsize=label_font_size)
        ax.tick_params(axis='y', which='both', labelsize=tick_font_size)
        ax.set_xlim(0, n_iters-1)
        plt.ticklabel_format(useOffset=False, axis='both')

    # Save the figure
    plt.savefig(fname, bbox_inches='tight', dpi=150)
    plt.savefig(fname[:-3] + 'png', bbox_inches='tight', dpi=150)
    plt.close()



def create_plot_posterior(params, plabs, cbars='red', nb=50, num=[]):


    fname = outdir+'/'+star+'_posterior.pdf'
    print('Creating ', fname)

    if (len(num) < 2):
        n = range(0, len(params))
    else:
        n = num

    priorf = prior_flags
    priorl = prior_vals

    plt.figure(1, figsize=(12, 4*(len(n))/n_columns_posterior))
    gs = gridspec.GridSpec(nrows=int(
        (len(n)+n_columns_posterior-1)/n_columns_posterior), ncols=n_columns_posterior)
    gs.update(wspace=0.025)
    j = int(0)
    for i in n:
        ax0 = plt.subplot(gs[j],rasterized=is_rasterized)
        vpar, lpar, rpar = find_vals_perc(params[i], 1.0)
        moda = my_mode(params[i])
        plt.axvline(x=vpar, c='#cc0000',label='Median',zorder=2,linewidth=1)
        plt.axvline(x=moda, c='#cc0000', ls='-.',label='Mode',zorder=2,linewidth=1)
        plt.xlabel(plabs[i])
        if (j % n_columns_posterior == 0):
            plt.ylabel('Frequency')
        plt.tick_params(axis='y', which='both',
                        direction='in', labelleft=False)
        plt.tick_params(axis='x', which='both', direction='in')
        counts, bins, patches = plt.hist(params[i], density=True, bins=nb, color=posterior_color,
                                 histtype='step', label='Posterior', alpha=1, linewidth=3)
        # Define the region to fill (from vpar-lpar to vpar+rpar)
        x_fill = np.linspace(vpar - lpar, vpar + rpar, 50)
        # Interpolate the histogram values to get the y-values for fill
        y_fill = np.interp(x_fill, (bins[:-1] + bins[1:]) / 2, counts)
        # Shade the area under the curve within the specified region
        plt.fill_between(x_fill, 0, y_fill, color=posterior_color, alpha=0.5,lw=0,label='68.3% C.I.')
        # Let us plot the prior ranges over the posterior distributions
        if is_plot_prior:
            lx, rx = ax0.get_xlim()
            if (priorf[i] == 'm' and lx < 0):
                lx = 1e-20  # avoid negative values
            locx = np.arange(lx, rx, (rx-lx)/1000.)
            lp = [None]*len(locx)
            for k in range(0, len(locx)):
                lp[k] = pti.get_priors(
                    priorf[i], [priorl[i*2], priorl[i*2+1]], locx[k])
            plt.plot(locx, lp, alpha=1, color=prior_color, label='Prior',linewidth=3)
        #
        plt.yticks([])
        #
        if (i == n[0]):
            plt.legend(loc=0, ncol=1, scatterpoints=1,
                       numpoints=1, frameon=True, fontsize=fos*0.5)
        j = int(j + 1)

    plt.savefig(fname, format='pdf', bbox_inches='tight',dpi=250)
    fname = outdir+'/'+star+'_posterior.png'
    plt.savefig(fname, format='png', bbox_inches='tight', dpi=300)
    plt.close()


def create_plot_correlation(params, plabs, col='red', mark='.', num=[],is_plot_prior=True,priorf=prior_flags,priorl=prior_vals):

    fname = outdir+'/'+star+'_correlations.pdf'
    print('Creating ', fname)

    if plot_kde_correlations:
        print("You set plot_kde_correlations=True, this may take a time to plot for runs with a lot of parameters")

    if (len(num) < 1):
        n = list(range(len(params)))
    else:
        n = num

    #Let us find the limits for each column of the plots
    limits = []
    for i in n:
        limits.append((min(params[i]),max(params[i])))


    plt.figure(1, figsize=(len(n), len(n)))
    nrows = len(n)
    ncols = len(n)
    gs = gridspec.GridSpec(nrows=nrows, ncols=ncols)
    gs.update(hspace=0.05,wspace=0.05)
    for o,i in enumerate(n):
        for p,j in enumerate(n):
            if j > i:
                continue
            ax0 = plt.subplot(gs[o*ncols+p],rasterized=is_rasterized)
            plt.tick_params(axis='y', which='both',
                            direction='in', labelleft=False)
            plt.tick_params(axis='x', which='both',
                            direction='in', labelbottom=False)
            plt.ticklabel_format(useOffset=False, axis='both')
            if (p == 0 and o > 0):
                plt.ylabel(plabs[i], fontsize=8)
                plt.tick_params(axis='y', which='both',
                                direction='in', labelleft=True,rotation=45,labelsize=6)
                plt.ticklabel_format(style='plain', axis='both', useOffset=False)
            if (i == n[len(n)-1]):
                plt.xlabel(plabs[j], fontsize=8)
                plt.tick_params(axis='x', which='both',
                                direction='in', labelbottom=True,rotation=45,labelsize=6)
                plt.ticklabel_format(style='plain', axis='both', useOffset=False)
            #PLOT POSTERIORS
            if j == i:
                vpar, lpar, rpar = find_vals_perc(params[i], 1.0)
                moda = my_mode(params[i])
                plt.axvline(x=vpar, c='#cc0000',label='Median',zorder=2,linewidth=1)
                plt.xlim(*limits[o])
                counts, bins, patches = plt.hist(params[i], density=True, bins=50, color=posterior_color,
                                 histtype='step', label='Posterior', alpha=1, linewidth=1)
                # Define the region to fill (from vpar-lpar to vpar+rpar)
                x_fill = np.linspace(vpar - lpar, vpar + rpar, 1000)
                # Interpolate the histogram values to get the y-values for fill
                y_fill = np.interp(x_fill, (bins[:-1] + bins[1:]) / 2, counts)
                # Shade the area under the curve within the specified region
                plt.fill_between(x_fill, 0, y_fill, color=posterior_color, alpha=0.5,lw=0,label='68.3% C.I.')
                # Let us plot the prior ranges over the posterior distributions
                if is_plot_prior:
                    lx, rx = ax0.get_xlim()
                    if (priorf[i] == 'm' and lx < 0):
                        lx = 1e-20  # avoid negative values
                    locx = np.arange(lx, rx, (rx-lx)/1000.)
                    lp = [None]*len(locx)
                for k in range(0, len(locx)):
                    lp[k] = pti.get_priors(
                            priorf[i], [priorl[i*2], priorl[i*2+1]], locx[k])
                plt.plot(locx, lp, alpha=1, color=prior_color, label='Prior',lw=1,zorder=2)
                plt.tick_params(axis='y',which='both', left=False, right=False)
                if j == n[0]: plt.legend(loc='upper right',bbox_to_anchor=(3.2, 0.92),fontsize=7)
            else:
                if plot_kde_correlations:
                    rindex = np.random.random_integers(0,len(params[j])-1,10000)
                    sns.kdeplot(params[j][rindex], params[i][rindex],levels=4,color='k')
                    plt.plot(params[j][rindex],params[i][rindex],'.',alpha=0.05,markersize=0.5,color='#006341')
                else:
                    z, xbins, ybins = np.histogram2d(params[j], params[i], bins=10)
                    plt.contourf(z.transpose(),extent=[xbins.min(),xbins.max(),ybins.min(),ybins.max()],cmap=correlation_cmap)
                plt.xlim(*limits[p])
                plt.ylim(*limits[o])

    plt.savefig(fname, format='pdf', bbox_inches='tight',dpi=250)
    fname = outdir+'/'+star+'_correlations.png'
    plt.savefig(fname, format='png', bbox_inches='tight', dpi=300)
    plt.close()


def create_corner_plot(params,labels,maxloglike,plot_parameters):
    import corner

    # update plot_parameters vector
    npp = list(plot_parameters)
    for o in range(len(plot_parameters)):
        npp[o] =  plot_parameters[o]

    # Let us take only the values to be plotted
    newpars = [0.0]*len(npp)
    newlabs = [0.0]*len(npp)
    true_params = [0.0]*len(npp)
    for o in range(0, len(npp)):
        newpars[o] = params[npp[o]]
        newlabs[o] = labels[plot_parameters[o]]
        true_params[o] = best_value(newpars[o], maxloglike, get_value)

    # Let us prepare the vector for corner
    data = np.zeros(shape=(len(newpars[0]), len(npp)))
    for o in range(len(newpars[0])):
        dumvec = []
        for m in range(len(npp)):
            dumvec.append(newpars[m][o])
        data[o] = dumvec

    figure = corner.corner(data, labels=newlabs,
                           quantiles=[0.16, 0.5, 0.84],
                           show_titles=True,)
    fname = outdir+'/'+star+'_corner.pdf'
    print('Creating ', fname)
    plt.savefig(fname, format='pdf', bbox_inches='tight',dpi=250)
    fname = outdir+'/'+star+'_corner.png'
    plt.savefig(fname, format='png', bbox_inches='tight', dpi=300)
    plt.close()
