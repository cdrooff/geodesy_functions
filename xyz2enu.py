def xyz2enu(refllh,xyz,covxyz=None):
    '''ECEF XYZ Coordinate system transformation to East-North-Up (m)

    Input: 1x3 reference lon-lat height vector at which the transformation matrix will be comupted
           Nx3 array of ECEF cartesian X-Y-Z
           Option to input a 3Nx3 X-Y-Z covariance matrix

    Output: Nx3 array of East-North-Up coordinates (m)
            Will return a 3Nx3 ENU covariance matrix if the XYZ covariance matrix is supplied as an argument

    '''

    import numpy as np
    import math

    # Define constants
    deg2rad = math.pi/180

    # Check to see if a covariance array has been passed in and set a flag
    if covxyz is None:
        do_covar = 0

    else:
        do_covar = 1


    if xyz.ndim == 1:
        nvec = 1
        enu = np.zeros(3)

        clat = math.cos(refllh[1]*deg2rad)
        slat = math.sin(refllh[1]*deg2rad)
        clon = math.cos(refllh[0]*deg2rad)
        slon = math.sin(refllh[0]*deg2rad)

        T = np.array(((-slon,clon,0),
                      (-slat*clon,-slat*slon,clat),
                      (clat*clon,clat*slon,slat)))

        enu = np.reshape(np.dot(T,np.reshape(xyz,(3,1))),(1,3))

        if do_covar:
            covenu = np.zeros((3,3))
            covenu = np.dot(np.dot(T,covxyz),np.transpose(T))


    else:

        nvec = np.shape(xyz)[0]

        if do_covar:
            covenu = np.zeros((3*nvec,3*nvec))

        enu = np.zeros((nvec,3))


        # Compute the transformation matrix
        clat = math.cos(refllh[1]*deg2rad)
        slat = math.sin(refllh[1]*deg2rad)
        clon = math.cos(refllh[0]*deg2rad)
        slon = math.sin(refllh[0]*deg2rad)

        T = np.array(((-slon,clon,0),
                      (-slat*clon,-slat*slon,clat),
                      (clat*clon,clat*slon,slat)))

        for ii in range(0,nvec):

            enu[ii,:] = np.reshape(np.dot(T,np.reshape(xyz[ii,:],(3,1))),(1,3))

        if do_covar:
            k2 = 3*ii
            covenu[ii:ii+3,k2:k2+3] = np.dot(np.dot(T,covxyz[ii:ii+3,k2:k2+3]),np.transpose(T))


    if do_covar:
        return enu, covenu

    else:
        return enu
