def filter_outliers(timeseries,sigma):
    '''
    Filters observations with high sigma values from a timeseries object

    Returns a timeseries object with observations removed for which the sigma value is exceeded on any component'''

    import copy
    import numpy as np

    out_timeseries = copy.deepcopy(timeseries)

    # Grab the indices for elements
    e_idx = list(np.where(timeseries.esig > sigma)[0])
    n_idx = list(np.where(timeseries.nsig > sigma)[0])
    h_idx = list(np.where(timeseries.hsig > sigma)[0])

    del_idx = e_idx

    for i in range(0,len(n_idx)):

        if n_idx[i] not in del_idx:
            del_idx.append(n_idx[i])

        elif h_idx[i] not in del_idx:
            del_idx.append(h_idx[i])

    # Sort the indices for deletion
    del_idx.sort()

    out_timeseries.time = np.delete(out_timeseries.time,del_idx)
    out_timeseries.llh = np.delete(out_timeseries.llh,del_idx,axis=0)
    out_timeseries.xyz = np.delete(out_timeseries.xyz,del_idx,axis=0)
    out_timeseries.esig = np.delete(out_timeseries.esig,del_idx)
    out_timeseries.nsig = np.delete(out_timeseries.nsig,del_idx)
    out_timeseries.hsig = np.delete(out_timeseries.hsig,del_idx)
    out_timeseries.enucov = np.delete(out_timeseries.enucov,del_idx,axis=0)


    # Re-compute east-north-up coordinates

    xyzdiff = out_timeseries.xyz - out_timeseries.xyz[0,:]
    refllh = out_timeseries.llh[0,:]
    enu_abs = xyz2enu(refllh,xyzdiff)
    out_timeseries.enu = enu_abs - enu_abs[0,:]


    out_timeseries.east = out_timeseries.enu[:,0]
    out_timeseries.north = out_timeseries.enu[:,1]
    out_timeseries.height = out_timeseries.enu[:,2]


    return out_timeseries
