#!/usr/bin/python
import numpy as num
import interpolator as interp
from interpolator.Interface import interpolateGenerator
import sys
t_offset = 0
filtstr = 'r'
mean_mag = 21.0

tms = []
ts = interpolateGenerator(sys.argv[1], mag0 = mean_mag, t0 = t_offset);
ts.getParams()
isper = ts.params['isPeriodic']
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
lc = interp.LightCurve(tms, isper)
lcn = lc.Realize(ra, dec, filtstr=filtstr, doAddErr = True, doDith = False, version="OpSim3_61")
#Available versions are: OpSim3_61, opsim1_29, opsim5_72, and cronos92
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
min = 50540
max = 50545
times = num.linspace(min,max,1000)
#dt = int(max - min)
#times = range(dt*10)
#times = num.array(times)/10. + min
mags = tsorig.evaluate(times, filt=filtstr)
print "# Values for the input light curve"
for i in range(len(times)):
    print times[i], mags[i]
print "\n"
