def window_timeseries(timeseries,t1,t2):
    '''Cuts a timeseries object to only lie within certain dates of interest

    Dates must be passed in float point or integer decimal year format.

    Entered dates do not need to precisely be dates of observations'''

    import copy
    import numpy as np

    ts_out = copy.deepcopy(timeseries)

    # Code expects t2 > t1
    # reassign variables if not fed in properly
    if (t1 > t2):
        t1,t2 = t2,t1

        # Grab indices of t1 and t2
    i_0 = int(min(np.where(timeseries.time >= t1)[0]))
    i_fin = int(max(np.where(timeseries.time <= t2)[0]))

    # Slice out the desired time window and return it as a new timeseries object
    ts_out.time = ts_out.time[i_0:i_fin]
    ts_out.llh = ts_out.llh[i_0:i_fin,:]
    ts_out.xyz =  ts_out.xyz[i_0:i_fin,:]

    ts_out.east = ts_out.east[i_0:i_fin]
    ts_out.north = ts_out.north[i_0:i_fin]
    ts_out.height = ts_out.height[i_0:i_fin]
    ts_out.enu = ts_out.enu[i_0:i_fin,:]

    # correct the east-north-up reference frame
    ts_out.east -= ts_out.east[0]
    ts_out.north -= ts_out.north[0]
    ts_out.height -= ts_out.height[0]
    ts_out.enu -= ts_out.enu[0,:]

    if hasattr(ts_out,'esig'):
        ts_out.esig = ts_out.esig[i_0:i_fin]

    if hasattr(ts_out,'nsig'):
        ts_out.nsig = ts_out.nsig[i_0:i_fin]

    if hasattr(ts_out,'hsig'):
        ts_out.hsig = ts_out.hsig[i_0:i_fin]

    if hasattr(ts_out,'enucov'):
        ts_out.enucov = ts_out.enucov[i_0:i_fin,:,:]

    return ts_out
