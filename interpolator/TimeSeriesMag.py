''' Create TimeSeriesMag object to hold time series infromation and meta data.
    Constructor takes an array of time points time, magnitude values values, 
    magnitude error valerr, a filter descriptor filterstr, a boolean which
    determines whether to calculate a spline for the time series calcspline,
    and a period for the time series in days period.  If the period argument
    is negative, the time series is assumed to be non periodic
    Modified:
    Jan 2008 by K. Simon Krughoff krughoff@astro.washington.edu
    March 2011 by K. Simon Krughoff for TVS group
'''
from Interpolate import makespline, evalper
from MagUtils import MagUtils
import numpy as num
class TimeSeriesMag:
    def __init__(self, time, values, valerrp, valerrm, filterstr, calcspline=True, period=-1.0, offset = 0, ra=None, dec=None, m5=None):
        ''' Construct TimeSeriesMag object from time sampling time, magnitude values values
            magnitude error values valerr, filter designation filterstr, whether to calculate
            a spline for the time series, and period of the time series in days.  If the 
            period argument is negative, the time series is assumed to be non periodic
            Inputs:
            time -- list of time values in days
            values -- list of magnitude values in magnitudes
            valerrp -- list of positive side photometric error values associated with each time point
            valerrm -- list of negative side photometric error values associated with each time point
            filterstr -- filter of time series (u, g, r, i, z, y)
            calcspline -- boolean, True to calculate a spline fit to the time series, False no spline calculated
            period -- period of the time series in days, should be negative if the time series is not periodic
            offset -- Offset in days of the beginning of this time series from MJD = 0
            ra -- RA in decimal degrees for this time series (default = None)
            dec -- DEC in decimal degrees for this time series (default = None)
            m5 -- If realized from a position in OpSim, this is the list of m5 values for each epoch
        '''
        #Set data for the object.
        self.getParams(time, values, valerrp, valerrm, m5, filterstr, calcspline, period, offset, ra, dec)

    def getParams(self, time, values, valerrp, valerrm, m5, filterstr, calcspline, period, offset, ra, dec):
        #Set time series arrays and meta data for a time series
        #Magnitude values
        values = num.asarray(values)
        #Time sampling for the time series
        time = num.asarray(time)
        #Magnitude errors
        valerrp = num.asarray(valerrp)
        valerrm = num.asarray(valerrm)
        self._magerrbright = valerrp
        self._magerrdim = valerrm
        self._mag = values
        self._m5s = m5
        self._time = time
        self._ra = ra
        self._dec = dec
        #Filter string accepted values are SDSS filter u, g, r, i, z, and y
        self._filter = filterstr
        #Set magnitude, magnitude error, time sampling in private arrays
        if calcspline:
            #Calculate a spline for the time series if calcspline is True
            self.setSpline(self._time, self.getFlux())
        else:
            #Set spline to None if calcspline is not True
            self._spline = None
        #Period of time series in days, negative if non periodic
        self._period = period
        #Offset of time series in days
        self._offset = offset
        #Check if magnitude and time arrays are consistent
        if(self._mag.any() and self._time.any()):
            assert len(self._mag) == len(self._time), "Magnitude array and time array must be the same length"

    def evaluate(self, epochs):
        return self.getSplineMags(epochs)

    def setSpline(self, time, flux):
        #Set spline for this time series
        self._spline = makespline(time, flux)

    def getSpline(self):
        #Get spline for this time series
        return self._spline

    def getTime(self):
        #Get time array for this time series
        return self._time 

    def getMag(self):
        #Get magnitude array for this time series
        return self._mag

    def getMagErrBright(self):
        #Get bright side magnitude error array for this time series
        return self._magerrbright

    def getMagErrDim(self):
        #Get dim side magnitude error array for this time series
        return self._magerrdim

    def getFlux(self):
        #Calculate and return the fluxes for the magnitude array in this time series.
        #Previous versions of this code have used Pogson magnitudes.  We have changed to 
        #ASINH mags since these behave well at negative fluxes.
        #See http://www.sdss.org/DR7/algorithms/fluxcal.html#counts2mag for a description and refs.
        fluxs = MagUtils().toFluxArr(self._mag, self._filter)
        return fluxs

    def getM5(self):
        return self._m5s

    def getFilter(self):
        #Get filter string for this time series
        return self._filter

    def getPeriod(self):
        #Get period of this time series
        return self._period
  
    def getSplineMags(self, dates):
        if(self._spline == None):
            print "Warning: Spline was not calculated on this time series"
            mags = None
        else:
            dates = num.asarray(dates)
            fluxs = num.asarray(evalper(self._spline, dates, self._offset, self._period))
            mags = MagUtils().toMagArr(fluxs, self._filter)
        return mags

    def getSplineFlux(self, dates):
        if(self._spline == None):
            print "Warning: Spline was not calculated on this time series"
            mags = None
        else:
            dates = num.asarray(dates)
            fluxs = num.asarray(evalper(self._spline, dates, self._offset, self._period))
        return fluxs 
