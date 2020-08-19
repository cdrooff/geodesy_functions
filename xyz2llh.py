def xyz2llh(xyz):
    '''ECEF XYZ coordinate Transformation to Lon-Lat-Height

    Input: Nx3 array of ECEF X-Y-Z Coordinates (m)

    Output: Nx3 array of lon-lat-height (deg, deg, m)'''

    import numpy as np
    import math

    # Define constants
    deg2rad = math.pi/180
    a = 6378137
    f = 1/298.257223563
    e2 = (2*f) - (f*f)

    if xyz.ndim == 1:
        nsites = 1

        llh = np.zeros(3)

        p = np.sqrt(np.dot(xyz[0:2],np.reshape(xyz[0:2],(2,1))))
        r = np.sqrt(np.dot(xyz[:],np.reshape(xyz[:],(3,1))))
        lon = math.atan2(xyz[1],xyz[0])

        # First iteration on lat and h
        # Assumes h=0

        lat = math.atan2(xyz[2]/p,(1-e2))
        n = a/np.sqrt((1-e2*math.sin(lat)**2))
        h = p/math.cos(lat) - n

        # Iterate until h converges

        oldh = -1e9
        itera = 0
        num = xyz[2]/p

        while ((abs(h-oldh) > 0.0001) and ((abs(lat) - math.pi/2) > 0.01)):
            itera = itera + 1
            oldh = h
            den = 1-e2*n/(n+h)
            lat = math.atan2(num,den)
            n = a/np.sqrt((1-e2*math.sin(lat)**2))
            h = p/math.cos(lat) - n

        llh[0] = lon/deg2rad + 360
        llh[1] = lat/deg2rad
        llh[2] = h


    else:
        nsites = np.shape(xyz)[0]


        # Initialize output array
        llh = np.zeros((nsites,3))


        for ii in range(0,nsites):
            p = np.sqrt(np.dot(xyz[ii,0:2],np.reshape(xyz[ii,0:2],(2,1))))
            r = np.sqrt(np.dot(xyz[ii,:],np.reshape(xyz[ii,:],(3,1))))
            lon = math.atan2(xyz[ii,1],xyz[ii,0])

            # First iteration on lat and h
            # Assumes h=0

            lat = math.atan2(xyz[ii,2]/p,(1-e2))
            n = a/np.sqrt((1-e2*math.sin(lat)**2))
            h = p/math.cos(lat) - n

            # Iterate until h converges

            oldh = -1e9
            itera = 0
            num = xyz[ii,2]/p

            while ((abs(h-oldh) > 0.0001) and ((abs(lat) - math.pi/2) > 0.01)):
                itera = itera + 1
                oldh = h
                den = 1-e2*n/(n+h)
                lat = math.atan2(num,den)
                n = a/np.sqrt((1-e2*math.sin(lat)**2))
                h = p/math.cos(lat) - n

            llh[ii,0] = lon/deg2rad
            llh[ii,1] = lat/deg2rad
            llh[ii,2] = h

    return llh
