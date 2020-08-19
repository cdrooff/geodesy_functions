def enu2xyz(refllh,enu):
    '''East-North-Up Coordinate System Transformation to ECEF XYZ

    Input: 1x3 reference lon-lat-height vector, Nx3 array of East-North-Up values transformed at refllh

    Output: Nx3 array of ECEF X-Y-Z coordinates

    '''

    import numpy as np
    import math

     # Define constants
    deg2rad = math.pi/180

    if enu.ndim == 1:
        xyz = np.zeros(3)

        clat = math.cos(refllh[1]*deg2rad)
        slat = math.sin(refllh[1]*deg2rad)
        clon = math.cos(refllh[0]*deg2rad)
        slon = math.sin(refllh[0]*deg2rad)

        T = np.array(((-slon,clon,0),
                      (-slat*clon,-slat*slon,clat),
                      (clat*clon,clat*slon,slat)))

        xyz = np.reshape(np.dot(T,np.reshape(enu,(3,1))),(1,3))

    else:

        nvec = np.shape(enu)[0]

        xyz = np.zeros((nvec,3))

        # Compute the transformation matrix
        clat = math.cos(refllh[1]*deg2rad)
        slat = math.sin(refllh[1]*deg2rad)
        clon = math.cos(refllh[0]*deg2rad)
        slon = math.sin(refllh[0]*deg2rad)

        T = np.array(((-slon,clon,0),
                      (-slat*clon,-slat*slon,clat),
                      (clat*clon,clat*slon,slat)))

        for ii in range(0,nvec):

            xyz[ii,:] = np.reshape(np.dot(T,np.reshape(enu[ii,:],(3,1))),(1,3))


    return llh2xyz(refllh) + xyz
