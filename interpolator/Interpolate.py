''' Wrapper class written by Andy Becker (UW) for the interpolation routines in scipy 
    makespline creates a spline with smoothing factor calculated from the dynamic range of the
    dependent variable
    evalper evaluates the input spline assuming a periodic input curve
    evalnonper evaluates the input spline assuming a non-periodic input curve
    splineinterp calls the appropriate evalutaion algorithm based on the boolean value isperiodic
    Modified:
    Jan 2008 by K. Simon Krughoff krughoff@astro.washington.edu
'''
#from scipy.interpolate import splrep, splev
from scipy.interpolate import UnivariateSpline
import numpy as num

def makespline(x, y, sfactor=0.01):
    ''' Make a spline interpolation of the input curve defined by x and y.
        Calculate smoothing factor, s, based on the dynamic range of the 
        dependent variable, y.
        Inputs:
        x -- independent variable for curve
        y -- dependent variable for curve
        sfactor -- smoothing factor
        Return:
        spline of input lightcurve
    '''
    val = None
    min = num.fabs(y).min()
    max = num.fabs(y).max()
    #The s value is used to constrain
    #the maximum allowable total square
    #difference between data and spline.
    #The sfactor is the average percent
    #deviation which I have set to 1%
    #(K. Krughoff Jan 2008)
    if y.min() == 0:
        val = (max - min)/2.
    else:
        val = min
    s = len(y)*(sfactor*val)**2
    tck = UnivariateSpline(x, y, s=s)
    return tck

def evalper(tck, xinterpolate, x0, xp):
    ''' Evaluate a periodic spline based on spline tck, and an array of independent values xinterpolate.
        Start may be shifted by offset x0.
        Period (in days) is given by xp.
        Inputs:
        tck -- spline to evaluate
        xinterpolate -- x locations for spine evaluation
        x0 -- offset in start location (in days for the lightcurve simulator)
        xp -- period in days of the lightcurve
        Return:
        list of interpolated flux values
    '''
    xp = float(xp)
    #Period in days 
    xnew = xinterpolate - x0
    #Apply offset to input independent values
    xnew %= xp
    #Find location in period at each x location
    xnew /= xp
    #Convert period location to fraction of period
    ynew = tck(xnew)
    #Interpolated y values
    return ynew

def evalnonper(tck, xspline, xinterpolate, x0, noflux=0.0):
    ''' Evaluate a non-periodic spline based on spline tck, 
        original x values for spline xspline, and an array of independent values xinterpolate.
        Start may be shifted by offset, x0.
        Value to return when there is no flux may be specified, noflux.
        Inputs:
        tck -- spline to evaluate
        xspline -- original x values for the spline
        xinterpolate -- x values for the evalutation for the spline (in days)
        x0 -- offset in days for interpolation of the spline
        noflux -- flux value to return if outside the bounds of the original spline
        Return:
        list of interpolated flux values
    '''
    xinterpolate = xinterpolate - x0
    #Apply offset to input independent values
    ynew = tck(xinterpolate)
    #Evaluate spilne at each x location in xinterpolate
    ynew = num.where(xinterpolate < xspline[0], noflux, ynew)
    ynew = num.where(xinterpolate > xspline[-1], noflux, ynew)
    #Return noflux value if x value is not in range of original spline range
    return ynew

def splineinterp(tck, xspline, xvals, x0, isperiodic, xp=-1):
    ''' Execute appropriate interpolation algorithm based on periodicity isperiodic
        Inputs:
        tck -- spline to evalutate
        xspline -- x values used to create the spline tck
        xvals -- x locations at which to evaluate the spline tck
        x0 -- offset in days to apply to the values xvals
        isperiodic -- boolean, True if the spline is periodic, False if not
        xp -- period of the input spline in days, should be negative if the spline is not periodic
        Return:
        list of interpolated flux values
    '''
    xspline = num.asarray(xspline)
    #numpy array of original x values for spline
    xvals = num.asarray(xvals)
    #numpy array of x values for interpolation
    if isperiodic:
        #execute evalper if curve is periodic
        return evalper(tck, xvals, x0, xp)
    else:
        #execte evalnonper if curve is not periodic
        return evalnonper(tck, xspline, xvals, x0, 0)
