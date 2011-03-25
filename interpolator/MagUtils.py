'''
    Created:
    Jan 2008 by K. Simon Krughoff krughoff@astro.washington.edu\
    Modified:
    March 2011 by K. Simon Krughoff fixes of docs.
'''
import numpy as num

class MagUtils:
    #This is just a guess at the softening parameter for the ASINH magnitude system.  
    #In reality, this should be related to the average sky and should vary with filter.
    #It corresponds to a zero flux mag of 27.5 and 1% deviation from Pogson at 25.
    _b = dict(u = 1.26e-10, g = 3.17e-11, r = 3.17e-11, i = 7.96e-11, z = 1.26e-10, y = 3.17e-10)
    #Zero magnitude flux density by definition is 3631Jy in the AB mag system
    _fo = 3631
    def getZeroPoint(self):
        #Get flux density zeropoint value for calculating fluxes for this time series
        return float(self._fo)
    def toFluxArr(self, mags, fs):
      b = None
      if(self._b.has_key(fs)):
        b = self._b[fs]
      else:
        b = 1.e-11
      return num.asarray(self._fo*2.*b*num.sinh(mags/(-2.5/num.log(10.)) - num.log(b)))
    def toMagArr(self, fluxs, fs):
      b = None
      if(self._b.has_key(fs)):
        b = self._b[fs]
      else:
        b = 1.e-11
      return num.asarray(-(2.5/num.log(10.))*(num.arcsinh((fluxs/self._fo)/(2.*b)) + num.log(b)))
    def toFlux(self, mag, fs):
      b = None
      if(self._b.has_key(fs)):
        b = self._b[fs]
      else:
        b = 1.e-11
      return self._fo*2.*b*num.sinh(mag/(-2.5/num.log(10.)) - num.log(b))
    def toMag(self, flux, fs):
      b = None
      if(self._b.has_key(fs)):
        b = self._b[fs]
      else:
        b = 1.e-11
      return -(2.5/num.log(10.))*(num.arcsinh((flux/self._fo)/(2.*b)) + num.log(b))

