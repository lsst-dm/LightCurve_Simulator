#!/usr/bin/python
import numpy as num
import interpolator as interp
#import pylab
import math
isper = 1
xp = 150 
t_offset = 0


infile = 'periodic.dat'
rfs = []
err = []
rms = []
time = []
for line in open(infile).readlines():

     t,rflux = map(float, line.split())


     rmag = -2.5*(math.log(rflux)/math.log(10)) + 35 
     err.append(1.0)
     rms.append(rmag)
     time.append(t)

tms = []
ts = interp.TimeSeriesMag(time, rms, err, err, 'g', period = xp);
tms.append(ts)
#save copy of the ideal light curve for comparison later
tsorig = ts

ra = []
dec = []

#put in the ra and dec pairs to sample
ra.append(210)
ra.append(210)
dec.append(-30)
dec.append(-20)

#send magnitudes and do the interpolation
lc = interp.LightCurve(isper, tms)
lcn = lc.Realize(ra, dec, t_offset, doAddErr = True, doDith = False, version="OpSim1_29")
#Available versions are: opsim1_29, opsim5_72, and cronos92
#If you wish to use the older version of the catalog, set the version to "Cronos92".
#Dithering can currently only be done on the older version.
#NOTE:  If resultant magerr = -9999 the m5 was brighter than the interpolated magnitude
timeres = []
for i in range(len(ra)):
    mytss = lcn[i].tss
    for ts in mytss:
        timeres = ts.getTime()
        magerrb = ts.getMagErrBright()
        magerrd = ts.getMagErrDim()
        magres = ts.getMag()
        m5res = ts.getM5()
        print "# Following are values for RA: %f, DEC: %f"%(ts._ra, ts._dec)
        for i in range(len(magres)):
            print timeres[i],magres[i],m5res[i],magerrb[i],magerrd[i]
        print "\n"

#Create a set of points to sample from the original light curve
min = min(timeres)
max = max(timeres)
dt = int(max - min)
times = range(dt)
times = num.array(times) + min
mags = tsorig.getSplineMags(times, t_offset)
fluxs = tsorig.getSplineFlux(times, t_offset)
print "# Values for the input light curve"
for i in range(len(times)):
    print times[i], mags[i]
print "\n"
