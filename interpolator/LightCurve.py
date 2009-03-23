''' Create a LightCurve object wich holds a set of TimeSeries objects and periodicity information
    The constructor takes a list of TimeSeries objects and a list of booleans which determine
    whether the TimeSeries objects are periodic.
    Realize takes a list of ra positions ra, dec positions dec, time offsets t_offset, whether to add
    random errors doAddErr, whether to use the dithered pointings from the database doDith
    Modified:
    March 2008 by K. Simon Krughoff krughoff@astro.washington.edu
'''
from Interpolate import splineinterp
from TimeSeriesMag import TimeSeriesMag
from MagUtils import *
from DB import DB
import numpy as num

class LightCurve:
    gamma = dict(u = 0.037, g = 0.038, r = 0.039, i = 0.039, z = 0.040, y = 0.040)
    #Constants to use in error calculation
    def __init__(self, isperiodic, ts):
        ''' 
        Construct LightCurve object with list of TimeSeries ts and list of booleans isperiodic 
        Inputs:
        isperiodic -- List of Booleans, True if the associated TimeSeriesMag object is periodic, False if not
        ts -- List of TimeSeriesMag objects for inclusion in the LightCurve object
        '''
        self.tss = ts
        self.isperiodic = isperiodic
    def Realize(self, ra, dec, t_offset, doAddErr = False, doDith = False, version="opsim1_29"):
        ''' 
        Realize the time sampling of a set of pointings based on a database of survey pointings 
        Inputs:
        ra -- list of RA values in degrees for sampling the lightcurves
        dec -- list of Declination values in degrees for sampling the lightcurves
        t_offset -- offset in days to apply to the time series when creating synthetic lightcurves
        doAddErr -- boolean, True to simulate errors on measurement based on 5 sigma limiting magnitude, False to return interpolated lightcurve uncorrected for errors
        doDith -- select sampling from the dithered version of the Cronos.92 pointings
        version -- select the version of the Operations Simulator to use.  Currently there are three options: Cronos92 (old), OpSim5_72, OpSim1_29
        Return:
        LightCurve object containing TimeSeries resampled based on the chosen operation simulator run
        '''
        ra = num.asarray(ra)
        #numpy array of ra positions
        dec = num.asarray(dec)
        #numpy array of dec positions
        toffset = []
        #time offsets
        scalartypes = (int, float, long)
        if isinstance(t_offset, scalartypes):
            #Determine if a list of offsets or simply a scalar offset were sent
            toffset = [t_offset]*len(ra)
        else:

            toffset = t_offset
        t_offset = num.asarray(toffset)
        #numpy array of time offsets
        assert len(ra) == len(dec), "ra and dec arrays must be the same length"
        assert len(t_offset) == len(ra), "t_offset must be a scalar or an array of the same length as ra and dec"
        #Check to make sure array lengths agree
        db = DB()
        lc = []
        for i in range(len(ra)):
            #Loop over positions for each TimeSeries
            tsi = []
            for ts in self.tss:
                tmpmag = []
                #Loop over TimeSeries array
                fs = ts.getFilter().lower()
                #Filter string
                if fs not in self.gamma:
                    #We only know sloan filters, so if something different assume r
                    print "Don't know parameters for filter",fs,"\n Assuming r..."
                    fs = 'r'
                (time, m5) = db.getTimeMagSQL(ra[i], dec[i], fs, doDith, version)
                #Get time sampling and 5 sigma limiting magnitude information from the database
                if len(time) == 0:
                    #If no data in database return a None object
                    tsi.append(TimeSeriesMag(None, None, None, None, fs, calcspline = False, ra = ra[i], dec = dec[i]))
                    continue
                fluxinterp = splineinterp(ts.getSpline(), ts.getTime(), time, t_offset[i], self.isperiodic, ts.getPeriod())
                fluxinterp = num.asarray(fluxinterp)
                #Interpolatd flux values based on time sampling from database
               
                #Calculate total photometric error from interpolated magnitudes and 5 sigma limiting magnitues.
                #Systematic error is assumed to be 0.01 magnitudes 
                m1flux = MagUtils().toFluxArr(m5, fs)/5.
                sysmag = (MagUtils().toMagArr(fluxinterp, fs) - 0.01)
                sysfluxerr = MagUtils().toFluxArr(sysmag,fs) - fluxinterp 
                sigs = []
                magerrfp = []
                magerrfm = []
                tmpflux = []
                if doAddErr:
                    #Add random error to interpolated magnitudes
                    tmpflux = fluxinterp + m1flux*num.random.normal(0,1,len(m1flux))
                    tmpflux = tmpflux + sysfluxerr*num.random.normal(0,1,len(sysfluxerr))
                    tmpmag = MagUtils().toMagArr(tmpflux, fs)
                    sigs = tmpflux/m1flux
                    
                else:
                    tmpflux = fluxinterp
                    tmpmag = MagUtils().toMagArr(tmpflux, fs)
                    sigs = tmpflux/m1flux

                x = 10.0**(0.4*(tmpmag-m5))
                magerr = num.sqrt((0.04 - self.gamma[fs])*x + self.gamma[fs]*x**2)
                magerr = num.sqrt(magerr**2 + 0.01**2)

                magerrfp = tmpmag - MagUtils().toMagArr(tmpflux + m1flux, fs)
                magerrfm = MagUtils().toMagArr(tmpflux - m1flux, fs) - tmpmag
                #at S/N > 10 the errors from flux differ from those calculated from the m5 by less
                #than the systematic error of 0.01 mag
                magerrfp = num.where(sigs > 10, magerr, magerrfp)
                magerrfm = num.where(sigs > 10, magerr, magerrfm)
                #if the detection is < 1 sigma indicate this by setting error to -9999
                magerrfp = num.where(sigs < 1, -9999, magerrfp)
                magerrfm = num.where(sigs < 1, -9999, magerrfm)
                tsi.append(TimeSeriesMag(time, tmpmag, magerrfp, magerrfm, fs, calcspline = False, ra = ra[i], dec = dec[i], m5 = m5))
                #Append resampled TimeSeries to output list of TimeSeries
            lc.append(LightCurve(self.isperiodic, tsi))
            #Append LightCurve for each ra/dec location 
        db.close()
        return lc
