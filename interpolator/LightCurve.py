''' Create a LightCurve object wich holds a set of TimeSeries objects and periodicity information
    The constructor takes a list of TimeSeries objects and a list of booleans which determine
    whether the TimeSeries objects are periodic.
    Realize takes a list of ra positions ra, dec positions dec, whether to add
    random errors doAddErr, whether to use the dithered pointings from the database doDith
    Modified:
    March 2008 by K. Simon Krughoff krughoff@astro.washington.edu
    March 2011 by K. Simon Krughoff for TVS group
'''
from Interpolate import splineinterp
from TimeSeriesMag import TimeSeriesMag
from MagUtils import MagUtils
from DB import DB
import numpy as num

class LightCurve:
    #Constants to use in error calculation
    gamma = dict(u = 0.037, g = 0.038, r = 0.039, i = 0.039, z = 0.040, y = 0.040)

    def __init__(self, ts, isperiodic):
        ''' 
        Construct LightCurve object with list of TimeSeries ts and list of booleans isperiodic 
        Inputs:
        ts -- List of TimeSeriesMag objects for inclusion in the LightCurve object
        isperiodic -- List of Booleans, True if the associated TimeSeriesMag object is periodic, False if not
        '''
        self.tss = ts
        self.isperiodic = isperiodic

    def Realize(self, ras, decs, doAddErr = False, doDith = False, version="opsim3_61"):
        ''' 
        Realize the time sampling of a set of pointings based on a database of survey pointings 
        Inputs:
        ras -- list of RA values in degrees for sampling the lightcurves
        decs -- list of Declination values in degrees for sampling the lightcurves
        doAddErr -- boolean, True to simulate errors on measurement based on 5 sigma limiting magnitude, 
                             False to return interpolated lightcurve uncorrected for errors
        doDith -- select sampling from the dithered version of the OpSim pointings
        version -- select the version of the Operations Simulator to use.  
                   Currently there are four options: Cronos92 (old), OpSim5_72, OpSim1_29, OpSim3_61 (default)
        Return:
        LightCurve object containing TimeSeries resampled based on the chosen operation simulator run
        '''
        #numpy array of ra positions
        ras = num.asarray(ras)
        #numpy array of dec positions
        decs = num.asarray(decs)
        assert len(ras) == len(decs), "ra and dec arrays must be the same length"
        db = DB()
        mutils = MagUtils()
        lc = []
        #Loop over positions for each TimeSeries
        for ra,dec in zip(ras,decs):
            tsi = []
            #Loop over TimeSeries array
            for ts in self.tss:
                tmpmag = []
                #Filter string
                fs = ts.getFilter().lower()
                if fs not in self.gamma:
                    #We only know sloan filters, so if something different assume r
                    print "Don't know parameters for filter",fs,"\n Assuming r..."
                    fs = 'r'
                #Get time sampling and 5 sigma limiting magnitude information from the database
                (time, m5) = db.getTimeMagSQL(ra, dec, fs, doDith, version)
                #If no data in database return a None object
                if len(time) == 0:
                    tsi.append(TimeSeriesMag(None, None, None, None, fs, calcspline = False, ra = ra, dec = dec))
                    continue
                #Interpolatd flux values based on time sampling from database
                fluxinterp = mutils.toFluxArr(ts.evaluate(time), fs)

                #Calculate total photometric error from interpolated magnitudes and 5 sigma limiting magnitues.
                #Systematic error is assumed to be 0.01 magnitudes 
                m1flux = mutils.toFluxArr(m5, fs)/5.
                sysmag = (mutils.toMagArr(fluxinterp, fs) - 0.01)
                sysfluxerr = mutils.toFluxArr(sysmag, fs) - fluxinterp 
                sigs = []
                magerrfp = []
                magerrfm = []
                tmpflux = []
                if doAddErr:
                    #Add random error to interpolated magnitudes
                    tmpflux = fluxinterp + m1flux*num.random.normal(0,1,len(m1flux))
                    tmpflux = tmpflux + sysfluxerr*num.random.normal(0,1,len(sysfluxerr))
                    tmpmag = mutils.toMagArr(tmpflux, fs)
                    sigs = tmpflux/m1flux
                    
                else:
                    tmpflux = fluxinterp
                    tmpmag = mutils.toMagArr(tmpflux, fs)
                    sigs = tmpflux/m1flux

                x = 10.0**(0.4*(tmpmag-m5))
                magerr = num.sqrt((0.04 - self.gamma[fs])*x + self.gamma[fs]*x**2)
                magerr = num.sqrt(magerr**2 + 0.01**2)

                magerrfp = tmpmag - mutils.toMagArr(tmpflux + m1flux, fs)
                magerrfm = mutils.toMagArr(tmpflux - m1flux, fs) - tmpmag

                #at S/N > 10 the errors from flux differ from those calculated from the m5 by less
                #than the systematic error of 0.01 mag
                magerrfp = num.where(sigs > 10, magerr, magerrfp)
                magerrfm = num.where(sigs > 10, magerr, magerrfm)

                #if the detection is < 1 sigma indicate this by setting error to -9999
                magerrfp = num.where(sigs < 1, -9999, magerrfp)
                magerrfm = num.where(sigs < 1, -9999, magerrfm)
                #Append resampled TimeSeries to output list of TimeSeries
                tsi.append(TimeSeriesMag(time, tmpmag, magerrfp, magerrfm, fs, calcspline = False, ra = ra, dec = dec, m5 = m5))
            #Append LightCurve for each ra/dec location 
            lc.append(LightCurve(tsi, self.isperiodic))
        db.close()
        return lc
