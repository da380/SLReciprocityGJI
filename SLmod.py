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
b    = 6371000.          # mean radius of the Earth
g    = 9.81              # mean surface gravity
G    = 6.674e-11         # gravitational constant
rhow = 1000.0            # density of water
rhoi =  917.0            # density of ice
CC = 8.038e37            # polar moment of inertia 
AA = 8.012e37            # equatorial moment of inertia
Om = 2.*pi/(24.*3600.)   # mean rotation rate
Me = 5.974e24            # mass of the Earth
##########################################################


##########################################################
# set some options
ep = 1.e-6               # tolerance for iterations
kind = 'complex'         # type of fields
csphase = -1             # Condon Shortley Phase
normalization = 'ortho'  # normalisation convention
##########################################################


######################################################################
# function to plot geographic data on GL grid. If complex data is
# input, only the real part is plotted
def geo_plot(data,cstring='RdBu'):
    ax = plt.axes(projection=ccrs.PlateCarree())
    plt.pcolormesh(data.lons()-180,data.lats(),np.real(data.data),shading = 'gouraud',cmap=cstring)
    ax.coastlines()
    plt.colorbar()
    plt.show()
    return


######################################################################
# function to read in present sea level and ice thickness and
# interpolate onto GL-grid (from ice6g0.dat)
def get_sl_ice_data(L):

    # initialise the grids
    sl = pysh.SHGrid.from_zeros(lmax=L,grid='GLQ',kind=kind)
    ice = sl.copy()
    
    dat = np.loadtxt('data/ice6g0.dat')

    nlat = int(dat[0,0])
    nlon = int(dat[0,1])
    
    sld = np.zeros([nlat,nlon],dtype = complex)
    latd = np.zeros([nlat,nlon])
    lond = np.zeros([nlat,nlon])
    icd = np.zeros([nlat,nlon],dtype = complex)

        
    k = 0
    for ilon in range(nlon):
        for ilat in range(nlat):
            k = k + 1
            latd[ilat,ilon] = dat[k,1]
            lond[ilat,ilon] = dat[k,0]+180.
            icd[ilat,ilon] = dat[k,2]
            sld[ilat,ilon] = -dat[k,3]
             
    
    # form the interpolating function for SL
    fun = interpolate.interp2d(lond[0,:],latd[:,0],np.real(sld),kind='linear')
    nlatgrid = sl.nlat
    nlongrid = sl.nlon
    lat = sl.lats()
    lon = sl.lons()
    
    # interpolate onto sl grid
    for ilat in range(nlatgrid):
        for ilon in range(nlongrid):
            sl.data[ilat,ilon] = fun(lon[ilon],lat[ilat])


    # form the interpolating function for ice
    fun = interpolate.interp2d(lond[0,:],latd[:,0],np.real(icd),kind='linear')
    

    #interpolate onto ice grid
    for ilat in range(nlatgrid):
        for ilon in range(nlongrid):
            ice.data[ilat,ilon] = fun(lon[ilon],lat[ilat])
            
    # correct sl where there is land-based ice
    #sl = sl + ice
    for ilat in range(nlatgrid):
        for ilon in range(nlongrid):
            if(ice.data[ilat,ilon] > 0 and sl.data[ilat,ilon] < 0):
                sl.data[ilat,ilon] = sl.data[ilat,ilon] + ice.data[ilat,ilon]
        
    return sl,ice               



###################################################################
# integrates a function over the surface
def surface_integral(fun):
    fun_lm = fun.expand(normalization=normalization,csphase=csphase)
    int = np.sqrt(4*pi)*b*b*np.real(fun_lm.coeffs[0,0,0])
    return int


###################################################################
# integrates a function over the oceans
def ocean_integral(C,fun):
    tmp    = C*fun
    tmp_lm = tmp.expand(normalization=normalization,csphase=csphase)
    int = np.sqrt(4*pi)*b*b*np.real(tmp_lm.coeffs[0,0,0])
    return int


###################################################################
# function that returns the ocean function and ocean area
def ocean_function(sl0,ice0):
    C = sl0.copy()
    nlat = sl0.nlat
    nlon = sl0.nlon
    
    # set the ocean function
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
def ocean_mask(C,val = np.nan):
    mask = C.copy()
    mask.data = np.where(np.real(C.data) == 1.,1.,val)
    return mask

#############################################################
# returns function equal to 1 on land and val in oceans
def land_mask(C,val = np.nan):
    mask = C.copy()
    mask.data = np.where(np.real(C.data) == 0.,1.,val)
    return mask


#############################################################
# returns function equal to 1 where there is ice and val
# elsewhere
def ice_mask(C,ice0,val = np.nan):
    mask = (1-C)*ice0
    mask.data = np.where(np.real(mask.data) > 0.,1.,val)
    return mask


#############################################################
# returns function equal to 1 where there is ice and val
# elsewhere
def greenland_mask(C,ice0,val = np.nan):
    mask = (1-C)*ice0
    mask.data = np.where(np.real(mask.data) > 0.,1.,val)
    nlat = mask.nlat
    nlon = mask.nlon    
    for ilat in range(nlat):
        lat = mask.lats()[ilat]
        for ilon in range(nlon):
            lon = mask.lons()[ilon]
            if(lat < 60 or lon < 50):
                mask.data[ilat,ilon] = val                
    return mask


#############################################################
# returns function equal to 1 where there is ice and val
# elsewhere
def antarctica_mask(C,ice0,val = np.nan):
    mask = (1-C)*ice0
    mask.data = np.where(np.real(mask.data) > 0.,1.,val)
    nlat = mask.nlat
    nlon = mask.nlon    
    for ilat in range(nlat):
        lat = mask.lats()[ilat]
        if(lat > 0):
            mask.data[ilat,:] = val                
    return mask


#############################################################
# function to read in the generalised Love numbers from a file
# and return them up to the required degree
def get_love_numbers(L):
    data = np.loadtxt('data/love_U.dat')
    hu = data[0:L+1,1]
    ku = data[0:L+1,3]
    data = np.loadtxt('data/love_P.dat')
    hp = data[0:L+1,1]
    kp = data[0:L+1,3]
    return hu,ku,hp,kp


#############################################################
# sets a spatial grid's values to zero in the southern
# hemisphere
def zero_southern_hemisphere(fin):
    nlat = fin.nlat
    lats = fin.lats()
    fout = fin.copy()
    for ilat in range(nlat):
        lat = lats[ilat]
        if(lat < 0.):
            fout.data[ilat,:] = 0.    
    return fout


#############################################################
# sets a spatial grid's values to zero in the northern
# hemisphere
def zero_northern_hemisphere(fin):
    nlat = fin.nlat
    lats = fin.lats()
    fout = fin.copy()
    for ilat in range(nlat):
        lat = lats[ilat]
        if(lat > 0.):
            fout.data[ilat,:] = 0.    
    return fout


#####################################################################
# returns the jx = i_{zx} and jy = i_{zy} components of the inertia
# tensor perturbation from the gravitational potential
def inertia_tensor_perturbation(phi_lm):
    j = np.zeros(2)
    j[0] = -np.sqrt(5./(6.*pi))*(b**3/G)*np.real(phi_lm.coeffs[0,2,1]) 
    j[1] = -np.sqrt(5./(6.*pi))*(b**3/G)*np.imag(phi_lm.coeffs[0,2,1])
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
def centrifugal_perturbation(om,psi_lm):
    psi_lm.coeffs[1,2,1] = b**2*Om*np.sqrt((2.*pi)/15.)*( om[0] + 1j*om[1])
    psi_lm.coeffs[0,2,1] = b**2*Om*np.sqrt((2.*pi)/15.)*(-om[0] + 1j*om[1])
    return


    
#####################################################################
# function to solve the fingerprint problem for a given direct load
def fingerprint(C,zeta,rotation=True):

    # get the maximum degree
    L = C.lmax

    # get the ocean area
    A = surface_integral(C)
    
    # get the love numbers
    data = np.loadtxt('data/love.dat')
    h = data[:L+1,1] + data[:L+1,3]
    k = data[:L+1,2] + data[:L+1,4]
    ht = data[2,5]
    kt = data[2,6]

    
    # calculate the average change in sea level
    slu = -surface_integral(zeta)/(rhow*A)
    onegrid = pysh.SHGrid.from_zeros(lmax=L,grid='GLQ',kind=kind)
    onegrid.data[:,:] = 1.
    sl = slu*onegrid
    
    # initialise displacement and potential perturbations
    u_lm   = pysh.SHCoeffs.from_zeros(lmax=sl.lmax,normalization=normalization, \
                                      kind=kind,csphase=csphase)
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
        sigma_lm  = sigma.expand(normalization=normalization,csphase=csphase)
        
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
            centrifugal_perturbation(om,psi_lm)

        # get the spatial fields
        u   =   u_lm.expand(grid='GLQ')
        phi = phi_lm.expand(grid='GLQ')
        psi = psi_lm.expand(grid='GLQ')


        # update the sea level
        fac = ocean_integral(C,g*u + phi + psi)/(g*A) + slu
        sl    = -1.*(g*u + phi + psi)/g + fac*onegrid

        it = it+1
        if(it == 0):
            slc = np.max(np.abs(sl.data))
        else:
            err = np.max(np.abs(sl.data - sl0.data))/np.abs(slc)
            print('iteration = ',it,'relative change = ',err)

        
        # store the most recent solution
        sl0 = sl.copy()
        
    if(not rotation):
        om = np.zeros(2)
        
    return sl,u,phi,om,psi




#####################################################################
def point_load(L,lats,lons,angle = 0.,w=[sentinel]):
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
               
    pl_lm = pysh.SHCoeffs.from_zeros(lmax=L,csphase=-1,normalization='ortho',kind='complex')
               
    for isource in range(len(lats)):

        ylm = pysh.expand.spharm(pl_lm.lmax,ths[isource],phs[isource], \
                                 csphase=csphase,normalization=normalization,kind=kind)
        
        for l in range(0,pl_lm.lmax+1):
            fac = np.exp(-l*(l+1)*t)
            pl_lm.coeffs[0,l,0] = pl_lm.coeffs[0,l,0] + w[isource]*ylm[0,l,0]*fac
            for m in range(1,l+1):
                pl_lm.coeffs[0,l,m] = pl_lm.coeffs[0,l,m] + (-1)**m*w[isource]*ylm[1,l,m]*fac
                pl_lm.coeffs[1,l,m] = pl_lm.coeffs[1,l,m] + (-1)**m*w[isource]*ylm[0,l,m]*fac


    pl_lm = (1/b**2)*pl_lm
    pl = pl_lm.expand(grid='GLQ')
               
    return pl





    
