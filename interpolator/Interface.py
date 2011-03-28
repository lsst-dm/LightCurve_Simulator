import sys
import numpy
import MySQLdb
from scipy.interpolate import InterpolatedUnivariateSpline

FILTERS = ['u', 'g', 'r', 'i', 'z', 'y']

#############
# BASE CLASS

class variabilityGenerator:
    def __init__(self):
        # attributes set during evaluate; may be float or numpy arrays
        self.reset()
        
    def reset(self):
        self.dMag = {}
        for filt in FILTERS:
            self.dMag[filt] = None

    def getParams(self):
        # must be overridden for every derived class
        pass
        
    def evaluate(self, epochs):
        # must be overridden for every derived class
        self.reset()
        pass

#############
# DERIVED CLASSES
#
# FOR INTERFACE WITH A TEXT FILE THAT WILL BE USED TO SPLINE-INTERPOLATE A LIGHTCURVE

class interpolateGenerator(variabilityGenerator):
    def __init__(self, infile, mag0 = 0., t0 = 0.):
        variabilityGenerator.__init__(self)
        self.filename = infile
        self.params   = {}
        self._t0      = t0
        self._m0      = mag0

    def getParams(self):
        # Required header format:
        # -File stars on the line below-
        # Period        : "0.358943659988 days"
        # isPeriodic    : True or False
        # Normalization : <g_mag>=0
        # Spectrum      : "Flat"; or filename
        # 
        # Phase u_mag  g_mag  r_mag  i_mag  z_mag  y_mag
        # XXX.X XXX.X  ...
        
        self.params['filename'] = self.filename
        fh = open(self.filename)
        period = fh.readline().split()[3]
        self.params['isPeriodic'] = fh.readline().split()[3]
        fh.close() 
        if self.params['isPeriodic']:
            self.params['period'] = float(period)
        else:
            self.params['lifetime'] = float(period)
        self.params['tOff']     = self._t0
        self.params['magOff']   = self._m0

        
    def evaluate(self, epochs, filt = None, isPerfect = True):
        lc    = numpy.loadtxt(self.params['filename'], unpack = True, comments='#')
        xval  = lc[0]
        uLc   = lc[1]+self.params['magOff']
        gLc   = lc[2]+self.params['magOff']
        rLc   = lc[3]+self.params['magOff']
        iLc   = lc[4]+self.params['magOff']
        zLc   = lc[5]+self.params['magOff']
        yLc   = lc[6]+self.params['magOff']

        splines = {}
        if isPerfect:
            # InterpolatedUnivariateSpline explicitly goes through each data point
            splines['u'] = InterpolatedUnivariateSpline(xval, uLc)
            splines['g'] = InterpolatedUnivariateSpline(xval, gLc)
            splines['r'] = InterpolatedUnivariateSpline(xval, rLc)
            splines['i'] = InterpolatedUnivariateSpline(xval, iLc)
            splines['z'] = InterpolatedUnivariateSpline(xval, zLc)
            splines['y'] = InterpolatedUnivariateSpline(xval, yLc)
        else:
            # UnivariateSpline smooths
            splines['u'] = UnivariateSpline(xval, uLc)
            splines['g'] = UnivariateSpline(xval, gLc)
            splines['r'] = UnivariateSpline(xval, rLc)
            splines['i'] = UnivariateSpline(xval, iLc)
            splines['z'] = UnivariateSpline(xval, zLc)
            splines['y'] = UnivariateSpline(xval, yLc)
            
        epochs -= self.params['tOff']
        
        if self.params.has_key('lifetime'):  
            spleval = epochs
        elif self.params.has_key('period'):
            spleval = epochs / self.params['period'] - epochs // self.params['period']
        else:
            raise Exception("No lifetime or period specified for this light curve")

        if filt == None:
            for f in FILTERS:
                self.dMag[f] = splines[f](spleval)
            return self.dMag
        else:
            assert(filt in FILTERS)
            self.dMag[filt] = splines[filt](spleval)
            return self.dMag[filt]

            

if __name__ == '__main__':
    # example of how to use an interpolator
    infile    = sys.argv[1]
    generator = interpolateGenerator(infile)
    generator.getParams()

    mjds      = numpy.arange(52000, 52010, 0.33)
    generator.evaluate(mjds)
    for i in xrange(len(mjds)):
        print ",".join([str(generator.dMag[f][i]) for f in FILTERS])

