#!/usr/bin/env python
import interpolator as interp
import scipy as num
import math
isper = 0

dmod = 20 # 20th mag at peak
t_offset = 52000 # realize with time offset of 52000 days to put it in the survey

infile = 'sn1a_lc.v1.1.dat'
Bfs = []
Vfs = []
err = []
Bms = []
Vms = []
time = []
for line in open(infile).readlines():

     t,U,B,V,R,I,J,H,K = map(float, line.split())


     err.append(1.0)
     Bms.append(B + dmod)
     Vms.append(V + dmod)
     time.append(t)

tms = []
#list of TimeSeries sending mags
ts = interp.TimeSeriesMag(time, Bms, err, err, 'g', offset=t_offset);
tms.append(ts)
ts = interp.TimeSeriesMag(time, Vms, err, err, 'r', offset=t_offset);
tms.append(ts)
tsorig = ts

ra = []
dec = []

#put in the ra and dec pairs to sample
ra.append(220)
ra.append(210)
dec.append(-30)
dec.append(60)

#send magnitudes and do the interpolation
lc = interp.LightCurve(tms, isper)
#Available versions are: opsim1_29, opsim5_72, and cronos92
#If you wish to use the older version of the catalog, set the version to "Cronos92".
#Dithering can currently only be done on the older version.
#NOTE:  If resultant mag = -66 it was outside the original lightcurve
#NOTE:  If resultant magerr = -9999 the m5 was brighter than the interpolated magnitude
#in this case, the m5 is returned instead of the interpolated magnitude
lcn = lc.Realize(ra, dec, doAddErr = True, doDith = False, version="OpSim3_61")
timeres = []
for i in range(len(ra)):
    mytss = lcn[i].tss
    for ts in mytss:
        print "#RA and Dec for this realization are: ",ts._ra, ts._dec
        if not ts.getTime().any():
           print "#This one has no time series info, probably because it is outside the survey region."
        else:
            timeres = ts.getTime()
            magerrbright = ts.getMagErrBright()
            magerrdim = ts.getMagErrDim()
            magres = ts.getMag()
            for i in range(len(magres)):
                print timeres[i],magres[i],magerrbright[i],magerrdim[i]
            print "\n"

#Spline and Interpolate original for comparison
min = min(timeres)
max = max(timeres)
dt = int(max - min)
times = range(dt)
times = num.array(times) + min
time = num.array(time)
mags = tsorig.evaluate(times)
fluxs = tsorig.getSplineFlux(times)
print "#Values for the interpolated input light curve"
for i in range(len(times)):
    print times[i], mags[i]
