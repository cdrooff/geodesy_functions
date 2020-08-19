def plot_timeseries(timeseries,show_errorbars=False):
    '''Reads in a timeseries object created by read_timeseries and plots out the timeseries on all three components

    Inputs:
    timeseries.time       Nx1 list of dates (decimal year)
    timeseries.sitename   1-4 character string with site name
    timeseries.comment    Comment string
    timeseries.llh        Nx3 array of lon-lat-height values
    timeseries.east       Series of east values
    timeseries.esig       Series of east sigmas
    timeseries.north      Series of north values
    timeseries.nsig       Series of north sigmas
    timeseries.height     Series of height values
    timeseries.hsig       Series of height sigmas
    timeseries.enu        Array of stored e-n-u values (m)
    timeseries.enucov     3x3xN array of Covariance Matrices
    timeseries.xyz        Nx3 array of x-y-z (m) in a geocentric cartesian coordinate system
    '''

    import matplotlib.pyplot as plt

    ts = timeseries


    # Plot the timeseries
    fig, axs = plt.subplots(3,1,figsize=(20,17))
    fig.suptitle(ts.sitename,fontsize=80, y=0.95)
    fig.subplots_adjust(hspace = .25, wspace=.001)


    # Plot Position in cm
    scale = 100

    #set dotsize
    dotsize = 70

    axs[0].tick_params(labelsize=16)
    axs[1].tick_params(labelsize=16)
    axs[2].tick_params(labelsize=16)


        ### Plot the East component positions and associated uncertainties ###

    if show_errorbars is True:
        axs[0].errorbar(ts.time,ts.east*scale,color='blue',ecolor='black',yerr=ts.esig*scale,linewidth=0.5,capsize=5,markersize=1,fmt='o',zorder=0)
        axs[1].errorbar(ts.time,ts.north*scale,color='black',yerr=ts.nsig*scale,linewidth=0.5,capsize=5,markersize=1,fmt='o',zorder=0)
        axs[2].errorbar(ts.time,ts.height*scale,color='black',yerr=ts.hsig*scale,linewidth=0.5,capsize=5,markersize=1,fmt='o',zorder=0)




    axs[0].scatter(ts.time,ts.east*scale,color='blue',s=dotsize,zorder=1)
    axs[0].set_ylabel('Longitude (cm)',size=30)


    ### Plot the North component positions and associated uncertainties ###


    axs[1].scatter(ts.time,ts.north*scale,color='green',s=dotsize,zorder=1)
    axs[1].set_ylabel('Latitude (cm)',size=30)


    ### Plot the Vertical component positions and associated uncertainties ###

    axs[2].scatter(ts.time,ts.height*scale,color='red',s=dotsize,zorder=1)
    axs[2].set_ylabel('Height (cm)',size=30)

    plt.xlabel('Year',size=30)

    return fig, axs
