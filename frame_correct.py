def frame_correct(timeseries,ref_frame='noam'):
    '''Corrects timeseries data to be relative to a desired reference frame
    By default assumes that the reference frame to correct for is stable North America as defined by Argus et al. (2010) in the GEODVEL modoel
    Can optionally pass in a separate tectonic unit as defined in GEODVEL using the ref_frame argument

    Works on a preexisting timeseries object and spits out a modified version of the input timeseries

    '''

    import numpy as np

    # calculate enu and xyz velocities and covariance matrices
    venu,vxyz,enucov,xyzcov = calc_geodvel('noam',timeseries.xyz[0])

    t0 = timeseries.time[0]
    tnorm = timeseries.time - t0

    timeseries.east -= venu[0,0]*tnorm
    timeseries.north -= venu[0,1]*tnorm
    timeseries.height -= venu[0,2]*tnorm

    timeseries.enu = np.transpose(np.array((timeseries.east,timeseries.north,timeseries.height)))

    timeseries.xyz = enu2xyz(timeseries.llh[0],timeseries.enu)
    timeseries.llh = xyz2llh(timeseries.xyz)

    return timeseries
