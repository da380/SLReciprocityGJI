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
ep = 1.e-8             # tolerance for iterations
##########################################################


######################################################################
# function to plot geographic data on GL grid.
def plot(fun,cstring='RdBu',contour = False, ncont = 6,title = '',marker = [],clim = []):
    ax = plt.axes(projection=ccrs.PlateCarree())
    if(contour):
        plt.contourf(fun.lons()-180,fun.lats(),fun.data,cmap=cstring,levels = ncont)
    else:
        plt.pcolormesh(fun.lons()-180,fun.lats(),fun.data,shading = 'gouraud',cmap=cstring)
    ax.coastlines()
    plt.colorbar()
#    plt.title(title)
    if(len(marker) == 2):
        lat = marker[0]
        lon = marker[1]
        plt.plot([lon],[lat],marker='o', markersize=5, color="red")
    if(len(clim) == 2):
        plt.clim(clim)
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

########################################################
# returns thw value of a field at a given point via
# spherical harmonic expansion
def point_evaluation(fun,lat,lon):
    fun_lm = fun.expand()
    return pysh.expand.MakeGridPoint(fun_lm.coeffs,lat,lon+180)


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
    for ilat in range(C.nlat):
        for ilon in range(C.nlon):            
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
# returns function equal to 1 where over antarctica and val
# elsewhere
def antarctica_mask(sl0,ice0,val = np.nan):
    mask = sl0.copy()
    for ilat,lat in enumerate(mask.lats()):
        for ilon,lon in enumerate(mask.lons()):            
            icel = ice0.data[ilat,ilon]            
            if(icel > 0. and lat < 0.): 
                mask.data[ilat,ilon] = 1.
            else:
                mask.data[ilat,ilon] = val
    return mask


#############################################################
# returns function equal to 1 where over antarctica and val
# elsewhere
def greenland_mask(sl0,ice0,val = np.nan):
    mask = sl0.copy()
    for ilat,lat in enumerate(mask.lats()):
        for ilon,lon in enumerate(mask.lons()):            
            icel = ice0.data[ilat,ilon]            
            if(icel > 0. and lat > 0. and lon > 110): 
                mask.data[ilat,ilon] = 1.
            else:
                mask.data[ilat,ilon] = val
    return mask


# returns function equal to 1 in oceans for lat in [lat1,lat2]
# and equal to val elsewhere
def altimetry_mask(sl0,ice0,lat1 = -66., lat2 = 66.,val = np.nan):
    mask = sl0.copy()
    for ilat,lat in enumerate(mask.lats()):
        for ilon,lon in enumerate(mask.lons()):
            sll  = sl0.data[ilat,ilon] 
            icel = ice0.data[ilat,ilon]            
            if(rhow*sll - rhoi*icel >= 0. and icel == 0 and lat >= lat1 and lat <= lat2): 
                mask.data[ilat,ilon] = 1.
            else:
                mask.data[ilat,ilon] = val
    return mask


#############################################################
# sets a spatial grid's values to zero in the southern
# hemisphere
def zero_southern_hemisphere(fin):
    fout = fin.copy()
    for ilat,lat in enumerate(fout.lats()):
        if(lat < 0.):
            fout.data[ilat,:] = 0.
    return fout


#############################################################
# sets a spatial grid's values to zero in the northern
# hemisphere
def zero_northern_hemisphere(fin):
    fout = fin.copy()
    for ilat,lat in enumerate(fout.lats()):
        if(lat > 0.):
            fout.data[ilat,:] = 0.
    return fout


#####################################################################
# returns the jx = Om*i_{zx} and jy = Om*i_{zy} components of the
# torque perturbation associated with the gravitational potential
def inertia_tensor_perturbation(phi_lm):
    j = Om*np.sqrt(5./(12*pi))*(b**3/G)*phi_lm.coeffs[:,2,1] 
    return j
 

#####################################################################
# returns the rotation vector perturbations given those for the
# torque vector, j.
def rotation_vector_perturbation(j):
    om = j/(CC-AA)
    return om


######################################################################
# returns the centrifugal potential perturbation in spherical harmonic
# domain given the rotation vector perturbation
def centrifugal_perturbation_coefficients(om):
    psi_2m = np.zeros([2,3])
    psi_2m[:,1] = b**2*Om*np.sqrt((4*pi)/15.)*om
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
# returns a vector kk such that (kk,om) is equal to the
# centrifugal potential perturbation at the given location
def centrifugal_perturbation_vector(om,lat,lon):
    kk = np.zeros(2)
    th = (90-lat)*pi/180
    ph = lon*pi/180
    kk[0] = Om*b*b*np.cos(th)*np.sin(th)*np.cos(ph)
    kk[0] = Om*b*b*np.cos(th)*np.sin(th)*np.sin(ph)
    return kk

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
# function to calculate uniform sea level rise from driect load
def bathtub(C,zeta):
    L = C.lmax
    A = surface_integral(C)
    return -surface_integral(zeta)/(rhow*A)


#####################################################################
# returns the solid earth deformation associated with a given load
def loading_response(sigma):
    L = sigma.lmax
    h,k,ht,kt = love_numbers(L)
    u_lm   = sigma.expand(normalization = 'ortho')
    phi_lm = u_lm.copy()
    for l in range(L+1):
        u_lm.coeffs[:,l,:]   *= h[l]
        phi_lm.coeffs[:,l,:] *= k[l]
    u   = u_lm.expand(grid = 'GLQ')
    phi = phi_lm.expand(grid = 'GLQ')    
    return u,phi

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
            u_lm.coeffs[:,l,:]   =  h[l]*sigma_lm.coeffs[:,l,:]    
            phi_lm.coeffs[:,l,:] =  k[l]*sigma_lm.coeffs[:,l,:]    
        

        # add in the centrifugal contribution
        u_lm.coeffs[:,2,:]   +=  ht*psi_lm.coeffs[:,2,:]
        phi_lm.coeffs[:,2,:] +=  kt*psi_lm.coeffs[:,2,:]

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
def generalised_fingerprint(C,zeta,zeta_u,zeta_phi,kk,rotation=True):

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
        u_lm.coeffs[:,2,:]   += ht*psi_lm.coeffs[:,2,:]
        phi_lm.coeffs[:,2,:] +=  kt*psi_lm.coeffs[:,2,:]

        # get the centrifugal potential perturbation
        if(rotation):
            j  = inertia_tensor_perturbation(phi_lm) - kk
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
# returns a point load at a given geographic location with optional
# inverse Laplacian smoothing
def point_load(L,lats,lons,grid = 'GLQ',angle = 0.,w=[sentinel]):


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
            pl_lm.coeffs[0,l,0] +=  w[isource]*ylm[0,l,0]*fac
            for m in range(1,l+1):
                pl_lm.coeffs[0,l,m] += w[isource]*ylm[0,l,m]*fac
                pl_lm.coeffs[1,l,m] += w[isource]*ylm[1,l,m]*fac


    pl_lm = (1/b**2)*pl_lm
    pl = pl_lm.expand(grid = 'GLQ')
               
    return pl



#########################################################################
# returns the vector -\int_{\partial M} [\mathbf{x} \times (\Bom \times
# \mathbf{x})] \zeta_{\psi} \dd S needed to deal with psi measurements
def rotation_vector_from_zeta_psi(zeta_psi):
    zeta_psi_lm = zeta_psi.expand(normalization='ortho')
    kk = np.zeros(2)
    for i in range(2):
        om = np.zeros(2)
        om[i] = 1.
        phi_2m = centrifugal_perturbation_coefficients(om)
        kk[i] = np.sum(phi_2m[:,:3]*zeta_psi_lm.coeffs[:,2,:3])*b*b
    return kk


        
##########################################################################
# returns the adjoint loads for a sea level measurement at 
def sea_level_load(L,lat,lon,grid = 'GLQ',angle = 1.):
    lats = np.full((1),lat)
    lons = np.full((1),lon)
    zeta     = point_load(L,lats,lons,angle = angle,grid = grid)
    zeta_u   = pysh.SHGrid.from_zeros(lmax=L,grid = grid)
    zeta_phi = pysh.SHGrid.from_zeros(lmax=L,grid = grid)
    kk       = np.zeros(2)
    return zeta,zeta_u,zeta_phi,kk


##########################################################################
# returns the adjoint loads for a sea level measurement at 
def displacement_load(L,lat,lon,grid = 'GLQ',angle = 1.):
    lats = np.full((1),lat)
    lons = np.full((1),lon)
    zeta     =  pysh.SHGrid.from_zeros(lmax=L,grid = grid)
    zeta_u   = -1*point_load(L,lats,lons,angle = angle,grid = grid)
    zeta_phi =  pysh.SHGrid.from_zeros(lmax=L,grid = grid)
    kk       = np.zeros(2)
    return zeta,zeta_u,zeta_phi,kk


    
##########################################################################
# returns the adjoint loads for a measurement of the (l,m)th real spherical
# harmonic coefficient of the gravitational potential perturbation
def potential_coefficient_load(L,l,m,grid = 'GLQ',remove_psi = True):
    zeta   = pysh.SHGrid.from_zeros(lmax=L,grid = grid)
    zeta_u = pysh.SHGrid.from_zeros(lmax=L,grid = grid)
    zeta_phi_lm =  pysh.SHCoeffs.from_zeros(lmax=L,normalization = 'ortho')
    if(m >= 0):
        zeta_phi_lm.coeffs[0,l,m]  = -g/b**2
    else:
        zeta_phi_lm.coeffs[1,l,-m] = -g/b**2
    zeta_phi = zeta_phi_lm.expand(grid = grid)
    if(remove_psi):
        kk = -rotation_vector_from_zeta_psi(zeta_phi)
    else:
        kk = np.zeros(2)
    return zeta,zeta_u,zeta_phi,kk



############################################################################
# returns the adjoint loads corresponding to a sea surface height measurement
# at a given latitude and longitude
def sea_altimetry_load(sl0,ice0,lat1 = -66,lat2 = 66,grid = 'GLQ',remove_psi = True):
    L = sl0.lmax
    zeta = altimetry_mask(sl0,ice0,lat1,lat2,val = 0.0)
    A = surface_integral(zeta)
    zeta = zeta/A
    zeta_u   = -1*zeta.copy()
    zeta_phi = pysh.SHGrid.from_zeros(lmax = L,grid=grid)
    if(remove_psi):
        kk = -rotation_vector_from_zeta_psi(zeta)
    else:
        kk = np.zeros(2)
    return zeta,zeta_u,zeta_phi,kk


############################################################################
# returns the load average corresponding to the averaging function, w,
# following the method of Wahr et al. (1998)
def GRACE_average_measurement(phi,w,LT = 0):
    L = phi.lmax
    if(LT == 0):
        LT = L
    _,k,_,_ = love_numbers(L)
    sigma_lm = phi.expand(normalization = 'ortho')
    for l in range(L+1):
        if(l >= 2 and l <= LT):
            sigma_lm.coeffs[:,l,:] *= 1./k[l]
        else:
            sigma_lm.coeffs[:,l,:] = 0.            
    sigma = sigma_lm.expand(grid='GLQ')
    return surface_integral(w*sigma)


############################################################################
# returns the adjoint loads corresponding to a GRACE load average defined
# with respec to the input function w. The truncation degree for the
# average can be input as an optional parameter
def GRACE_average_load(w,LT = 0):
    L = w.lmax
    if(LT  == 0):
        LT = L
    zeta   = pysh.SHGrid.from_zeros(lmax = L,grid='GLQ')
    zeta_u = pysh.SHGrid.from_zeros(lmax = L,grid='GLQ')
    _,k,_,_ = love_numbers(L)
    w_lm = w.expand(normalization = 'ortho')    
    for l in range(L+1):
        if(l >= 2 and l <= LT):
            w_lm.coeffs[:,l,:] *= -g/k[l]
        else:
            w_lm.coeffs[:,l,:] = 0.            
    zeta_phi = w_lm.expand(grid='GLQ')
    kk = -rotation_vector_from_zeta_psi(zeta_phi)
    return zeta,zeta_u,zeta_phi,kk


####################################################
# returns the averaging function of Jekeli (1981)
# as discussed in Wahr et al (1998). Note that
# the averaging length, r, is input in km
def circular_averaging_function(L,r,lat0,lon0):
    th0 = (90-lat0)*pi/180
    ph0 = (lon0-180)*pi/180
    c = np.log(2)/(1-np.cos(1000*r/b))
    fac = 2*pi*(1-np.exp(-2*c))
    fac = c/(b*b*fac)
    w = pysh.SHGrid.from_zeros(lmax=L,grid = 'GLQ')
    for ilat,lat in enumerate(w.lats()):
        th = (90-lat)*pi/180
        fac1 = np.cos(th)*np.cos(th0)
        fac2 = np.sin(th)*np.sin(th0)
        for ilon,lon in enumerate(w.lons()):
            ph = lon*pi/180
            calpha = fac1 + fac2*np.cos(ph-ph0)
            w.data[ilat,ilon] = fac*np.exp(-c*(1-calpha))
    w_lm = w.expand()
    w_lm.coeffs[:,:2,:] = 0.
    w = w_lm.expand(grid = 'GLQ')
    return w



##############################################################
# sets Laplacian-type covariance
def laplace_covariance(L,std = 1.,s = 0.,mu = 0.):
    Q = np.zeros(L+1)
    norm = 0.
    for l in range(L+1):
        fac = 1.+mu*mu*l*(l+1)
        fac = fac**(-s)
        norm += (2*l+1)*fac/(4*pi)
        Q[l] = fac
    Q = b*b*std*std*Q/norm
    return Q

###############################################################
# generates random field with Laplacian like covariance
def random_field(Q):
    L = Q.size-1
    fun_lm = pysh.SHCoeffs.from_zeros(lmax=L,normalization = 'ortho')
    for l in range(L+1):
        fac = np.sqrt(Q[l])/b
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
        fac = Q[l]/(b*b)
        cf_lm.coeffs[0,l,0] = cf_lm.coeffs[0,l,0] + ylm[0,l,0]*fac
        for m in range(1,l+1):
            cf_lm.coeffs[0,l,m] += ylm[0,l,m]*fac
            cf_lm.coeffs[1,l,m] += ylm[1,l,m]*fac
    cf = cf_lm.expand(grid = 'GLQ')
    return cf


###################################################################
# returns a random ice model derived from a zero-mean rotationally
# invariant Gassian random field with the prescribed covariance

def random_ice_model(sl0,ice0,Q):
    return random_field(Q)*ice_mask(sl0,ice0,val = 0.)

###################################################################
# returns a random ocean model derived from a zero-mean rotationally
# invariant Gassian random field with the prescribed covariance
# note that the integral of the field over the oceans is set equal to
# zero in accordance with conservation of mass
def random_ocean_model(sl0,ice0,Q):
    C = ocean_function(sl0,ice0)
    rf = random_field(Q)*C
    rf -= C*surface_integral(rf)/surface_integral(C)
    return rf
    

####################################################
# returns the action of the covariance operator, Q,
# for a rotationally invariant Gaussian random field
# on an input function
def covariance_action(Q,fun):
    grid = fun.grid
    L = fun.lmax
    fun_lm = fun.expand(normalization = 'ortho')
    for l in range(L+1):
        fun_lm.coeffs[:,l,:] *= Q[l]
    return fun_lm.expand(grid = grid)




