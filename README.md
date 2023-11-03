
## Summary

Supporting material for [[1]](#1)
which has been accepted for publication in Geophysical Journal International. These codes are sufficient to reproduce all calculations and figures
within the  paper. This can be done either with the Jupyter notebook, MakeFigures.ipynb, or using the Python script, MakeFigures.py,
which saves the figures to a chosen directory. 

Functions needed to solve the (generalised) fingerprint problem are contained in SLmod.py. Also included is RFmod.py which contains functions
for generating Gaussian random fields on a sphere; these are, optionally, used within the main scripts to test the methods against
random changes of ice thickness.

Values for the present day sea level and ice thickness have been taken from the
model ICE6G of [[2]](#2) and [[3]](#3) . See https://www.atmosp.physics.utoronto.ca/~peltier/data.php for more details.

Tidal and (generalised) loading Love numbers have been pre-calculated for PREM [[4]](#4) using a Fortran code that can be found
within  https://github.com/da380/gia3D.git. 

If you make use of these codes within your own work, an acknowledgement would be appreciated along with a citation to [[1]](#1). 

## Dependencies

A list of all key import statements can be seen in the script SLmod.py. Beyond standard things (matplotlib, numpy, etc.), it is necessary to  install
the following packagees:

1. **pyshtools** which is used for fast spherical harmonic transformations. See https://shtools.github.io/SHTOOLS/using-with-python.html for details
on its use and installation.

2. **cartopy** which is used for plotting coastlines on various figures.  See https://scitools.org.uk/cartopy/docs/latest/ for details on its use and
installation.



## References
<a id="1">[1]</a> 
Al-Attar D., Syvret F., Crawford O., Mitrovica J.X., and Lloyd A.J., 2023.
Reciprocity and sensitivity kernels for sea level fingerprints. Submitted to *Geophys. J. Int.*.

<a id="2">[2]</a> 
Argus, D.F., Peltier, W.R., Drummond, R. and Moore, A.W.(2014) The Antarctica component of postglacial rebound model ICE-6G_C (VM5a) based upon GPS positioning, exposure age dating of ice thicknesses, and relative sea level histories. Geophys. J. Int., 198(1), 537-563, doi:10.1093/gji/ggu140.



<a id="3">[3]</a> 
Peltier, W.R., Argus, D.F. and Drummond, R. (2015) Space geodesy constrains ice-age terminal deglaciation: The global ICE-6G_C (VM5a) model. J. Geophys. Res. Solid Earth, 120, 450-487, doi:10.1002/2014JB011176.

<a id="4">[4]</a> 
Dziewonski, A.M. and Anderson, D.L., 1981. Preliminary reference Earth model. Physics of the earth and planetary interiors, 25(4), pp.297-356.
