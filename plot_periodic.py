#!/usr/bin/python
import numpy as num
import interpolator as interp
from interpolator.Interface import interpolateGenerator
import sys
import pylab

infile    = sys.argv[1]  # e.g. sampleRRly.txt
filter    = 'r'          # 'u', 'g', 'r', 'i', 'z', 'y'
t_offset  = 0.0          # days
mean_mags = [20, 21, 22, 23, 24, 25]

# where do you want to realize the lightcurves?
ra     = [210, 210]  # deg
decl   = [-30, -20]  # deg

# all the lightcurves you want sampled at each ra and decl
allModels = []
for mean_mag in mean_mags:
    mdata = interpolateGenerator(sys.argv[1], mag0 = mean_mag, t0 = t_offset) # read the input data
    mdata.getParams() # initialize it
    isper = mdata.params['isPeriodic']
    allModels.append(mdata)

# Available lcSim versions are: OpSim3_61, opsim1_29, opsim5_72, and cronos92
# Note all models have to be periodic or not, here
lcsim  = interp.LightCurve(allModels, isper)
allLcs = lcsim.Realize(ra, decl, filtstr=filter, doAddErr = True, doDith = False, version="OpSim3_61")

for nPos in range(len(ra)):
    lcs = allLcs[nPos]
    for nMag in range(len(lcs.tss)):
        mean_mag    = mean_mags[nMag]
        period      = allModels[nMag].params['period']
        realization = lcs.tss[nMag]

        # single values
        raLc     = realization._ra
        declLc   = realization._dec

        # arrays
        dates    = realization.getTime()
        mags     = realization.getMag()
        magerrs1 = realization.getMagErrBright()
        magerrs2 = realization.getMagErrDim()

        idx      = num.where((magerrs1 > 0) & (magerrs2 > 0))
        if len(idx[0]) == 0:
            continue

        dates    = dates[idx]
        mags     = mags[idx]
        magerrs1 = magerrs1[idx]
        magerrs2 = magerrs2[idx]

        fig = pylab.figure()
        sp1 = fig.add_subplot(211)
        sp1.errorbar(dates, mags, 0.5 * (magerrs1 + magerrs2), fmt = 'ro')
        sp1.set_xlabel('MJD')
        sp1.set_ylabel('Mag')
        sp1.set_title('Mean Mag = %f; Ra = %f; Decl = %f' % (mean_mag, raLc, declLc))
        ymin, ymax = sp1.get_ylim()
        sp1.set_ylim((ymax, ymin))

        sp2 = fig.add_subplot(212)
        phase = dates / period - dates // period 
        sp2.errorbar(phase, mags, 0.5 * (magerrs1 + magerrs2), fmt = 'ro')
        sp2.set_xlabel('Phase')
        sp2.set_ylabel('Mag')
        ymin, ymax = sp2.get_ylim()
        sp2.set_ylim((ymax, ymin))

        # theory curve; just place it in sp2 since it would be massive in sp1
        phases = num.arange(0, 1, 0.01)
        times  = period * phases
        mags   = allModels[nMag].evaluate(times, filt=filter)
        sp2.plot(phases, mags, 'k--')

pylab.show()

