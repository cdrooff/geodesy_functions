def calc_geodvel(inplate,xyz):
    '''
    Given an input position vector computes the ENU velocity based on GEODVEL
    (Argus et al., 2010, doi: 10.1111/j.1365-246X.2009.04463.x)
    A Geocenter correction for use with ITRF2008 updated from that paper is included

    Inputs:
    inplate - a 4 character lower case string specifying the major tectonic plate with which to compute the velocity of
    xyz - an Nx3 array of ECEF XYZ coordinates at which to compute plate velocity

    Outputs:
    Venu - a local 3xN East-North-Up velocity vector
    Vcov - 3x3xN covariance array in East-North-Up
    Tcov - 3x3xN covariance array in ECEF X-Y-Z
    '''

    import numpy as np

    # Load the GEODVEL model
    plates, omegas, omegacov, geocenter_xyz, geocenter_cov = load_geodvel_model()

    # Grab the index for the correct plate
    try:
        iplate = plates.index(inplate)
    except ValueError:
        print("""calc_geodvel error! Cannot load plate data. Plate name needs to be a 4-character lowercase string
        Valid entires are:
        'anta'
        'arab'
        'aust'
        'eura'
        'noam'
        'indi'
        'nazc'
        'nubi'
        'pcfc'
        'soam'
        'soma'
        'carb'""")


    om = omegas[:,iplate]
    k = 3*iplate
    omcov = omegacov[k:k+3,k:k+3]


    # Check if xyz is a vector
    if xyz.ndim == 1:
        nsites = 1
        rmat = np.zeros((3*nsites,3))

        rmat = ((0, xyz[2], -xyz[1]), (-xyz[2], 0, xyz[0]), (xyz[1], -xyz[0], 0))


    else:
        nsites = np.shape(xyz)[0]
        rmat = np.zeros((3*nsites,3))


        for i in range(0,nsites):
            k = 3*i

            rmat[k:k+3,:] = ((0, xyz[i,2], -xyz[i,1]), (-xyz[i,2], 0, xyz[i,0]), (xyz[i,1], -xyz[i,0], 0))


    vxyz = np.dot(rmat,om)
    tcov = np.dot(np.dot(rmat,omcov),np.transpose(rmat))

    ### Correct for Geocenter Motion ###
    gc_rmat = np.zeros((3*nsites,3))

    for i in range(0,nsites):
        k = 3*i
        gc_rmat[k:k+3,:] = np.eye(3)

    vxyz -= np.dot(gc_rmat,geocenter_xyz)
    tcov += np.dot(gc_rmat,np.dot(geocenter_cov,np.transpose(gc_rmat)))


    venu = np.zeros(np.shape(xyz))
    vcov = np.zeros((3,3*nsites))


    llh = xyz2llh(xyz)

    if llh.ndim == 1:
        venu, covenu = xyz2enu(llh,vxyz,tcov)

    else:
        for i in range(0,nsites):
            k = 3*i

            #print(vxyz[k:k+3])

            venu[i,:], covenu = xyz2enu(llh[i,:],vxyz[k:k+3],tcov[k:k+3,k:k+3])

            vcov[:,k:k+3] = covenu

        vxyz = np.reshape(vxyz,(nsites,3))

    return venu, vxyz, vcov, tcov
