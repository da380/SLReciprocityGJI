
# Summary

Supporting material for the paper "Reciprocity and sensitivity kernels for sea level fingerprints" by Al-Attar *et al.* (2023)
which has been submitted to GJI. The codes and model files provided here are sufficient to reproduce all calculations and figures
within the submitted paper. The Jupyter notebook MakeFigures.ipynb can be used to generate the figures interactively or to experiment more
generally with the methods developed within the paper. The Python3 script MakeFigures.py can be used to generate all figures and save them to
a given location.

# Dependencies

A list of all key import statements can be seen in the script SLmod.py. Beyond standard things (matplotlib, numpy), it is necessary to have installed
the following packagees:

1. **pyshtools** which is used for fast spherical harmonic transformations. See https://shtools.github.io/SHTOOLS/using-with-python.html for details
on its use and installation.

2. **cartopy** which is used for plotting coastlines on various figures.  See https://scitools.org.uk/cartopy/docs/latest/ for details on its use and
installation. 




