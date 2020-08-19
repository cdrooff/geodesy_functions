def read_timeseries(pfile): #, sigtol, refvel, refllh):
    '''Read in a pfile and output a timeseries object. Timeseries attributes listed below.

    timeseries.time       Nx1 list of dates (decimal year)
    timeseries.sitename   1-4 character string with site name
    timeseries.llh        Nx3 array of lon-lat-height values
    timeseries.xyz        Nx3 array of x-y-z (m) in an ECEF coordinate system
    timeseries.east       Series of east values
    timeseries.esig       Series of east sigmas
    timeseries.north      Series of north values
    timeseries.nsig       Series of north sigmas
    timeseries.height     Series of height values
    timeseries.hsig       Series of height sigmas
    timeseries.enu        Nx3 array of stored e-n-u values (m)
    timeseries.enucov     3x3xN array of Covariance Matrices
    '''

    import numpy as np

    class timeseries():

        def __init__(self):

            self.time = list(np.genfromtxt(pfile)[:,0])
            self.sitename = np.loadtxt(pfile, dtype=str)[0,1]
            self.llh = np.genfromtxt(pfile)[:,3:6]

            self.xyz = llh2xyz(self.llh)

            # Calculate East-North-Up (m) displacements at all timesteps
            xyzdiff = self.xyz - self.xyz[0,:]
            refllh = self.llh[0,:]
            enu_abs = xyz2enu(refllh,xyzdiff)
            self.enu = enu_abs - enu_abs[0,:]

            self.east = self.enu[:,0]
            self.north = self.enu[:,1]
            self.height = self.enu[:,2]

            self.esig = np.genfromtxt(pfile)[:,6]/1000
            self.nsig = np.genfromtxt(pfile)[:,7]/1000
            self.hsig = np.genfromtxt(pfile)[:,8]/1000

            corr_en = np.genfromtxt(pfile)[:,9]
            corr_ev = np.genfromtxt(pfile)[:,10]
            corr_nv = np.genfromtxt(pfile)[:,11]

            self.enucov = np.zeros((len(self.east), 3, 3))

            # Create an array which represent the covariance matrix for each set of observations
            for ii in range(0,len(corr_en)):

                # Initialize the covariance array
                cov = np.zeros((3,3))

                # Populate the diagonal elements
                self.enucov[ii,0,0] = self.esig[ii]**2
                self.enucov[ii,1,1] = self.nsig[ii]**2
                self.enucov[ii,2,2] = self.hsig[ii]**2

                # Populate the non-diagonal elements
                self.enucov[ii,0,1] = corr_en[ii]*self.esig[ii]*self.nsig[ii]
                self.enucov[ii,0,2] = corr_ev[ii]*self.esig[ii]*self.hsig[ii]
                self.enucov[ii,1,2] = corr_nv[ii]*self.nsig[ii]*self.hsig[ii]

                # Matrix is symmetric
                self.enucov[ii,1,0] = self.enucov[ii,0,1]
                self.enucov[ii,2,0] = self.enucov[ii,0,2]
                self.enucov[ii,2,1] = self.enucov[ii,1,2]


    tseries = timeseries()

    # Define a sigma tolerance for use in filtering timeseries
    sigtol = 0.04

    # Filter the output timeseries object
    ts_out = filter_outliers(tseries,sigtol)

    return ts_out
