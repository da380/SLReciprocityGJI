
# Summary

Supporting material for the paper "Reciprocity and sensitivity kernels for sea level fingerprints" by Al-Attar *et al.* (2023)
which has been submitted to GJI. The codes and model files provided are sufficient to reproduce all calculations and figures
within the submitted paper. The Jupyter notebook MakeFigures.ipynb can be used to generate the figures interactively or to experiment more
generally with the methods developed within the paper. The Python3 script MakeFigures.py can be used to generate all figures and save them to
a given location. The computational routines needed to solve the sea level equation are located in SLmod.py. There is also a module RFmod.py which
is used to generate Gaussian random fields on a sphere. This is, optionally, used within the main scripts to generate random ice thickness changes.

These codes have been written in collaboration between David Al-Attar and Frank Syvret.

# Dependencies

A list of all key import statements can be seen in the script SLmod.py. Beyond standard things (matplotlib, numpy, etc.), it is necessary to  install
the following packagees:

1. **pyshtools** which is used for fast spherical harmonic transformations. See https://shtools.github.io/SHTOOLS/using-with-python.html for details
on its use and installation.

2. **cartopy** which is used for plotting coastlines on various figures.  See https://scitools.org.uk/cartopy/docs/latest/ for details on its use and
installation.

# Data files

Within the data subdirectory there are three files. Two are the present-day sea level and ice thickness that have been taken from the
model ICE6G of Argus *et al.* (2014). See https://www.atmosp.physics.utoronto.ca/~peltier/data.php for the source of this data. The final
file contains generalised Love numbers for PREM of Dziewonski & Anderson (1981). The code used for these calculations can
be found within https://github.com/da380/gia3D.git. 




