import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import pyshtools as pysh
from numpy import pi as pi
from scipy import interpolate 
sentinel = object()


if __name__ == "__main__":
    pass


#########################################################
# set some constants
b    = 6368000.          # mean radius of the Earth
g    = 9.825652323       # mean surface gravity
G    = 6.6723e-11        # gravitational constant
rhow = 1000.0            # density of water
rhoi =  917.0            # density of ice
CC = 8.038e37            # polar moment of inertia 
AA = 8.012e37            # equatorial moment of inertia
Om = 2.*pi/(24.*3600.)   # mean rotation rate
Me = 5.974e24            # mass of the Earth
##########################################################


##########################################################
# set some options
ep = 1.e-6              # tolerance for iterations
##########################################################


######################################################################
# function to plot geographic data on GL grid. If complex data is
# input, only the real part is plotted
def plot(fun,cstring='RdBu'):
    ax = plt.axes(projection=ccrs.PlateCarree())
    plt.pcolormesh(fun.lons()-180,fun.lats(),fun.data,shading = 'gouraud',cmap=cstring)
    ax.coastlines()
    plt.colorbar()
    plt.show()
    return


######################################################################
# function to read in present sea level and ice thickness and
# interpolate onto GL-grid (from ice6g0.dat)
def get_sl_ice_data(L):

    # initialise the grids
    sl = pysh.SHGrid.from_zeros(lmax=L,grid = 'GLQ')
    ice = sl.copy()
    
    dat = np.loadtxt('data/ice6g0.dat')

    nlat = int(dat[0,0])
    nlon = int(dat[0,1])
    
    sld = np.zeros([nlat,nlon])
    latd = np.zeros([nlat,nlon])
    lond = np.zeros([nlat,nlon])
    icd = np.zeros([nlat,nlon])

        
    k = 0
    for ilon in range(nlon):
        for ilat in range(nlat):
            k = k + 1
            latd[ilat,ilon] =  dat[k,1]
            lond[ilat,ilon] =  dat[k,0]+180.
            icd[ilat,ilon]  =  dat[k,2]
            sld[ilat,ilon]  = -dat[k,3]
             
    
    # form the interpolating function for SL
    fun = interpolate.interp2d(lond[0,:],latd[:,0],sld,kind='linear')
    nlatgrid = sl.nlat
    nlongrid = sl.nlon
    lat = sl.lats()
    lon = sl.lons()
    
    # interpolate onto sl grid
    for ilat in range(nlatgrid):
        for ilon in range(nlongrid):
            sl.data[ilat,ilon] = fun(lon[ilon],lat[ilat])


    # form the interpolating function for ice
    fun = interpolate.interp2d(lond[0,:],latd[:,0],icd,kind='linear')
    

    #interpolate onto ice grid
    for ilat in range(nlatgrid):
        for ilon in range(nlongrid):
            ice.data[ilat,ilon] = fun(lon[ilon],lat[ilat])
            
    # correct sl where there is land-based ice
    for ilat in range(nlatgrid):
        for ilon in range(nlongrid):
            if(ice.data[ilat,ilon] > 0 and sl.data[ilat,ilon] < 0):
                sl.data[ilat,ilon] = sl.data[ilat,ilon] + ice.data[ilat,ilon]
        
    return sl,ice               



###################################################################
# integrates a function over the surface
def surface_integral(fun):
    fun_lm = fun.expand(lmax_calc = 0,normalization = 'ortho')
    int = np.sqrt(4*pi)*b*b*fun_lm.coeffs[0,0,0]
    return int


###################################################################
# integrates a function over the oceans
def ocean_integral(C,fun):
    tmp    = C*fun
    tmp_lm = tmp.expand(lmax_calc = 0,normalization = 'ortho')
    int = np.sqrt(4*pi)*b*b*tmp_lm.coeffs[0,0,0]
    return int


###################################################################
# function that returns the ocean function and ocean area
def ocean_function(sl0,ice0):
    C = sl0.copy()
    nlat = sl0.nlat
    nlon = sl0.nlon
    for ilat in range(nlat):
        for ilon in range(nlon):            
            sll  = sl0.data[ilat,ilon] 
            icel = ice0.data[ilat,ilon]            
            if(rhow*sll - rhoi*icel >= 0.): 
                C.data[ilat,ilon] = 1.
            else:
                C.data[ilat,ilon] = 0.
    return C


#############################################################
# returns function equal to 1 in oceans and val on land
def ocean_mask(sl0,ice0,val = np.nan):
    mask = sl0.copy()
    for ilat in range(mask.nlat):
        for ilon in range(mask.nlon):            
            sll  = sl0.data[ilat,ilon] 
            icel = ice0.data[ilat,ilon]            
            if(rhow*sll - rhoi*icel >= 0.): 
                mask.data[ilat,ilon] = 1.
            else:
                mask.data[ilat,ilon] = val
    return mask

#############################################################
# returns function equal to 1 on land and val in oceans
def land_mask(sl0,ice0,val = np.nan):
    mask = sl0.copy()
    for ilat in range(mask.nlat):
        for ilon in range(mask.nlon):            
            sll  = sl0.data[ilat,ilon] 
            icel = ice0.data[ilat,ilon]            
            if(rhow*sll - rhoi*icel >= 0.): 
                mask.data[ilat,ilon] = val
            else:
                mask.data[ilat,ilon] = 1.
    return mask


#############################################################
# returns function equal to 1 where there is ice and val
# elsewhere
def ice_mask(sl0,ice0,val = np.nan):
    mask = sl0.copy()
    for ilat in range(mask.nlat):
        for ilon in range(mask.nlon):            
            sll  = sl0.data[ilat,ilon] 
            icel = ice0.data[ilat,ilon]            
            if(rhow*sll - rhoi*icel < 0. and icel > 0.): 
                mask.data[ilat,ilon] = 1.
            else:
                mask.data[ilat,ilon] = val
    return mask


#############################################################
# sets a spatial grid's values to zero in the southern
# hemisphere
def zero_southern_hemisphere(fin):
    fout = fin.copy()
    ilat = 0
    for lat in fout.lats():
        if(lat < 0.):
            fout.data[ilat,:] = 0.
        ilat += 1
    return fout


#############################################################
# sets a spatial grid's values to zero in the northern
# hemisphere
def zero_northern_hemisphere(fin):
    fout = fin.copy()
    ilat = 0
    for lat in fout.lats():
        if(lat > 0.):
            fout.data[ilat,:] = 0.
        ilat += 1
    return fout


#####################################################################
# returns the jx = i_{zx} and jy = i_{zy} components of the inertia
# tensor perturbation from the gravitational potential
def inertia_tensor_perturbation(phi_lm):
    j = np.zeros(2)
    j[0] =  -np.sqrt(5./(12*pi))*(b**3/G)*phi_lm.coeffs[0,2,1] 
    j[1] =  -np.sqrt(5./(12*pi))*(b**3/G)*phi_lm.coeffs[1,2,1]
    return j


#####################################################################
# returns the rotation vector perturbations given those for the
# inertia tensor
def rotation_vector_perturbation(j):
    om = np.zeros(2)
    om = Om*j/(CC-AA)
    return om


######################################################################
# returns the centrifugal potential perturbation in spherical harmonic
# domain given the rotation vector perturbation
def centrifugal_perturbation_coefficients(om):
    psi_2m = np.zeros([2,3])
    psi_2m[0,1] = b**2*Om*np.sqrt((4*pi)/15.)*om[0]
    psi_2m[1,1] = b**2*Om*np.sqrt((4*pi)/15.)*om[1]
    return psi_2m


########################################################################
# returns the point value of the centrifugal potential perturbation at
# a given location
def centrifugal_perturbation_value(om,lat,lon):
    th = (90-lat)*pi/180
    ph = lon*pi/180
    psi = Om*b*b*np.cos(th)*np.sin(th)*(om[0]*np.cos(ph) + om[1]*np.sin(ph))
    return psi

########################################################################
# returns a vector lv such that (lv,om) is equal to the
# centrifugal potential perturbation at the given location
def centrifugal_perturbation_vector(om,lat,lon):
    lv = np.zeros(2)
    th = (90-lat)*pi/180
    ph = lon*pi/180
    lv[0] = Om*b*b*np.cos(th)*np.sin(th)*np.cos(ph)
    lv[0] = Om*b*b*np.cos(th)*np.sin(th)*np.sin(ph)
    return lv

################################################
# reads in and returns love numbers from a file
def love_numbers(L):
    data = np.loadtxt('data/love.dat')
    h = data[:L+1,1] + data[:L+1,3]
    k = data[:L+1,2] + data[:L+1,4]
    ht = data[2,5]
    kt = data[2,6]
    return h,k,ht,kt


################################################
# reads in and returns love numbers from a file
def generalised_love_numbers(L):
    data = np.loadtxt('data/love.dat')
    h_u   = data[:L+1,1]
    k_u   = data[:L+1,2]
    h_phi = data[:L+1,3]
    k_phi = data[:L+1,4]
    h     = h_u + h_phi
    k     = k_u + k_phi
    ht = data[2,5]
    kt = data[2,6]
    return h_u,k_u,h_phi,k_phi,h,k,ht,kt

    



#####################################################################
# function to solve the fingerprint problem for a given direct load
def fingerprint(C,zeta,rotation=True):

    # get the maximum degree
    L = C.lmax

    # get the ocean area
    A = surface_integral(C)
    
    # get the love numbers
    h,k,ht,kt = love_numbers(L)
    
    # calculate the average change in sea level
    slu = -surface_integral(zeta)/(rhow*A)
    onegrid = pysh.SHGrid.from_zeros(lmax=L,grid = 'GLQ')
    onegrid.data[:,:] = 1.
    sl = slu*onegrid
    
    # initialise displacement and potential perturbations
    u_lm   = pysh.SHCoeffs.from_zeros(lmax=sl.lmax,normalization = 'ortho')
    phi_lm = u_lm.copy()
    psi_lm = u_lm.copy()
        
    # store the initial guess
    sl0 = sl.copy()
    err = 1.

    # start the iterations
    it = -1
    while(err > ep):
        
        # compute the current loads
        sigma = rhow*C*sl + zeta
        sigma_lm  = sigma.expand(normalization = 'ortho')
        
        # determine the response to the loading
        for l in range(L+1):
            u_lm.coeffs[:,l,:]   =    h[l]*sigma_lm.coeffs[:,l,:]    
            phi_lm.coeffs[:,l,:] =    k[l]*sigma_lm.coeffs[:,l,:]    
        

        # add in the centrifugal contribution
        u_lm.coeffs[:,2,:]   =   u_lm.coeffs[:,2,:] + ht*psi_lm.coeffs[:,2,:]
        phi_lm.coeffs[:,2,:] = phi_lm.coeffs[:,2,:] + kt*psi_lm.coeffs[:,2,:]

        # get the centrifugal potential perturbation
        if(rotation):
            j  = inertia_tensor_perturbation(phi_lm)
            om = rotation_vector_perturbation(j)
            psi_2m = centrifugal_perturbation_coefficients(om)
            psi_lm.coeffs[:,2,:3] = psi_2m

        # get the spatial fields
        u   =   u_lm.expand(grid = 'GLQ')
        phi = phi_lm.expand(grid = 'GLQ')
        psi = psi_lm.expand(grid = 'GLQ')


        # update the sea level
        fac = ocean_integral(C,g*u + phi + psi)/(g*A) + slu
        sl    = -1.*(g*u + phi + psi)/g + fac*onegrid

        it = it+1
        if(it == 0):
            slnorm = np.max(np.abs(sl.data))
        else:
            err = np.max(np.abs(sl.data - sl0.data))/np.abs(slnorm)
            print('iteration = ',it,'relative change = ',err)

        
        # store the most recent solution
        sl0 = sl.copy()
        
    if(not rotation):
        om = np.zeros(2)
        
    return sl,u,phi,om,psi




#####################################################################
# function to solve the fingerprint problem for a given direct load
def generalised_fingerprint(C,zeta,zeta_u,zeta_phi,lv,rotation=True):

    # get the maximum degree
    L = C.lmax

    # get the ocean area
    A = surface_integral(C)
    
    # get the love numbers
    h_u,k_u,h_phi,k_phi,h,k,ht,kt = generalised_love_numbers(L)
        
    # calculate the average change in sea level
    slu = -surface_integral(zeta)/(rhow*A)
    onegrid = pysh.SHGrid.from_zeros(lmax=L,grid='GLQ')
    onegrid.data[:,:] = 1.
    sl = slu*onegrid
    
    # initialise displacement and potential perturbations
    u_lm   = pysh.SHCoeffs.from_zeros(lmax=sl.lmax,normalization = 'ortho')
    phi_lm = u_lm.copy()
    psi_lm = u_lm.copy()


    # expand the generalised loads
    zeta_u_lm   = zeta_u.expand(normalization = 'ortho')
    zeta_phi_lm = zeta_phi.expand(normalization = 'ortho')

    # start the iterations
    err = 1.
    it = -1
    while(err > ep):
        
        # compute the current loads
        sigma       = rhow*C*sl + zeta
        sigma_lm    = sigma.expand(normalization = 'ortho')
        
        # determine the response to the loading
        for l in range(L+1):
            u_lm.coeffs[:,l,:]   =    h[l]*sigma_lm.coeffs[:,l,:]        \
                                    + h_u[l]*zeta_u_lm.coeffs[:,l,:]     \
                                    + h_phi[l]*zeta_phi_lm.coeffs[:,l,:]
            phi_lm.coeffs[:,l,:] =    k[l]*sigma_lm.coeffs[:,l,:]        \
                                    + k_u[l]*zeta_u_lm.coeffs[:,l,:]     \
                                    + k_phi[l]*zeta_phi_lm.coeffs[:,l,:]


        # add in the centrifugal contribution
        u_lm.coeffs[:,2,:]   =   u_lm.coeffs[:,2,:] + ht*psi_lm.coeffs[:,2,:]
        phi_lm.coeffs[:,2,:] = phi_lm.coeffs[:,2,:] + kt*psi_lm.coeffs[:,2,:]

        # get the centrifugal potential perturbation
        if(rotation):
            j  = inertia_tensor_perturbation(phi_lm) + lv/Om
            om = rotation_vector_perturbation(j)
            psi_2m = centrifugal_perturbation_coefficients(om)
            psi_lm.coeffs[:,2,:3] = psi_2m
            
        # get the spatial fields
        u   =   u_lm.expand(grid = 'GLQ')
        phi = phi_lm.expand(grid = 'GLQ')
        psi = psi_lm.expand(grid = 'GLQ')

        # update the sea level
        fac = ocean_integral(C,g*u + phi + psi)/(g*A) + slu
        sl    = -1.*(g*u + phi + psi)/g + fac*onegrid

        it = it+1
        if(it == 0):
            slnorm = np.max(np.abs(sl.data))
        else:
            err = np.max(np.abs(sl.data - sl0.data))/np.abs(slnorm)
            print('iteration = ',it,'relative change = ',err)

        
        # store the most recent solution
        sl0 = sl.copy()
        
    if(not rotation):
        om = np.zeros(2)
        
    return sl,u,phi,om,psi





#####################################################################
def point_load(L,lats,lons,grid = 'GLQ',angle = 0.,w=[sentinel]):
# function returns a point load at a given geographic location
# Note that the result is defined on a unit sphere and so will
# need to be normalised in other cases

    if(len(lats)!=len(lons)):
        raise SystemExit('lats and lons are different sizes!')
        
    if(w[0]==sentinel):
        w = np.ones(len(lats))

    th = 0.4*angle*pi/180
    t  = th*th
        
    ths = 90.- lats
    phs = lons + 180.
               
    pl_lm = pysh.SHCoeffs.from_zeros(lmax=L,normalization = 'ortho')
               
    for isource in range(len(lats)):

        ylm = pysh.expand.spharm(pl_lm.lmax,ths[isource],phs[isource],normalization = 'ortho')
        
        for l in range(0,pl_lm.lmax+1):
            fac = np.exp(-l*(l+1)*t)
            pl_lm.coeffs[0,l,0] = pl_lm.coeffs[0,l,0] + w[isource]*ylm[0,l,0]*fac
            for m in range(1,l+1):
                pl_lm.coeffs[0,l,m] = pl_lm.coeffs[0,l,m] + w[isource]*ylm[0,l,m]*fac
                pl_lm.coeffs[1,l,m] = pl_lm.coeffs[1,l,m] + w[isource]*ylm[1,l,m]*fac


    pl_lm = (1/b**2)*pl_lm
    pl = pl_lm.expand(grid = 'GLQ')
               
    return pl


############################################################
# transforms a function under the Sobolev mapping
def sobolev_mapping(u,s,mu):   
    u_lm = u.expand(normalization = 'ortho')
    L = u_lm.lmax
    for l in range(L+1):
        fac = 1.0 + mu*mu*l*(l+1)
        fac = fac**(-s)
        u_lm.coeffs[:,l,:] = fac*u_lm.coeffs[:,l,:]
    u = u_lm.expand(grid = 'GLQ')
    return u
    

#############################################################
# computes the Sobolev inner product. Note that default
# values of s and mu given the special case of the L2 product
def inner_product(u,v,s = 0., mu = 0.):
    u_lm = u.expand(normalization = 'ortho')
    v_lm = u.expand(normalization = 'ortho')
    L = u_lm.lmax
    p = 0.
    for l in range(L+1):
        fac = 1.0+mu*mu*l*(l+1)
        fac = fac**s
        p = p + fac*u_lm.coeffs[0,l,0]*v_lm.coeffs[0,l,0]
        for m in range(1,l+1):
            p = p + fac*u_lm.coeffs[0,l,m]*v_lm.coeffs[0,l,m]
            p = p + fac*u_lm.coeffs[1,l,m]*v_lm.coeffs[1,l,m]            
    return p*b*b
        
        
##############################################################
# sets Laplacian-type covariance
def laplace_covariance(L,std = 1.,s = 0.,mu = 0.):
    Q = np.zeros(L+1)
    norm = 0.
    for l in range(L+1):
        fac = 1.+mu*mu*l*(l+1)
        fac = fac**(-s)
        norm = norm + (2*l+1)*fac/(4*pi)
        Q[l] = fac
    Q = std*std*Q/norm
    return Q

###############################################################
# generates random field with Laplacian like covariance
def random_field(Q):
    L = Q.size-1
    fun_lm = pysh.SHCoeffs.from_zeros(lmax=L,normalization = 'ortho')
    for l in range(L+1):
        fac = np.sqrt(Q[l])
        ranp = np.random.normal(size = l+1)
        rann = np.random.normal(size = l+1)
        fun_lm.coeffs[0,l,0] = fac*ranp[0]
        for m in range(1,l+1):
            fun_lm.coeffs[0,l,m] = fac*ranp[m]
            fun_lm.coeffs[1,l,m] = fac*rann[m]
    return fun_lm.expand(grid = 'GLQ')        

################################################################
# returns the two-point correlation function for a random field
# given the covariance operator
def correlation_function(Q,lat0 = 0.,lon0 = 0.):
    L = Q.size-1    
    ths = 90.- lat0
    phs = lon0 + 180.               
    cf_lm = pysh.SHCoeffs.from_zeros(lmax=L,normalization = 'ortho')               
    ylm = pysh.expand.spharm(cf_lm.lmax,ths,phs,normalization = 'ortho')        
    for l in range(L+1):
        fac = Q[l]
        cf_lm.coeffs[0,l,0] = cf_lm.coeffs[0,l,0] + ylm[0,l,0]*fac
        for m in range(1,l+1):
            cf_lm.coeffs[0,l,m] = cf_lm.coeffs[0,l,m] + ylm[0,l,m]*fac
            cf_lm.coeffs[1,l,m] = cf_lm.coeffs[1,l,m] + ylm[1,l,m]*fac
    cf = cf_lm.expand(grid = 'GLQ')
    return cf



