def llh2xyz(llh):
    ''' Lon-Lat-Height to ECEF X-Y-Z Coordinate System Transformation

    Input: Nx3 array of lon-lat-height and

    Output: an Nx3 array of x-y-z (m) in an ECEF reference frame'''

    import numpy as np
    import math

    # Define constants
    deg2rad = math.pi/180
    a = 6378137
    f = 1/298.257223563
    e2 = (2*f) - (f*f)

    if llh.ndim == 1:
        xyz = np.zeros(3)

        clat = math.cos(llh[1]*deg2rad)
        slat = math.sin(llh[1]*deg2rad)
        clon = math.cos(llh[0]*deg2rad)
        slon = math.sin(llh[0]*deg2rad)

        n = a/np.sqrt(1-e2*(slat**2))

        xyz[0] = (n + llh[2])*clat*clon
        xyz[1] = (n + llh[2])*clat*slon
        xyz[2] = (n*(1-e2) + llh[2])*slat

    else:
        nvec = np.shape(llh)[0]
        xyz = np.zeros((nvec,3))

        for ii in range(0,len(llh)):
            clat = math.cos(llh[ii,1]*deg2rad)
            slat = math.sin(llh[ii,1]*deg2rad)
            clon = math.cos(llh[ii,0]*deg2rad)
            slon = math.sin(llh[ii,0]*deg2rad)

            n = a/np.sqrt(1-e2*(slat**2))

            xyz[ii,0] = (n + llh[ii,2])*clat*clon
            xyz[ii,1] = (n + llh[ii,2])*clat*slon
            xyz[ii,2] = (n*(1-e2) + llh[ii,2])*slat

    return xyz
