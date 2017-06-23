import numpy as np
import glob
import matplotlib.pyplot as plt
from sedfit import source

for file in glob.glob('sedfit/M11_SEDs/Indiv/*.dat'):

    print file
    wav, flux = np.loadtxt(file, \
                            comments='#', \
                            unpack=True)
    o = np.where(wav < 500.)
    wav = wav[o]
    flux = flux[o]

    flux = flux/np.max(flux)

    src = source(0,0,0.1)
    src.sed.setWav(wav)
    src.sed.setPL(alpha=2.8, plnorm=0.2, turnover=40)
    src.sed.setBB(bbnorm=1, temp=40., beta=1.5)

    err = flux/5.
    o = np.where(wav < 35.)
    err[o] = 1.0*err[o]

    result = src.sed.fit(flux, err)

    result.params.valuesdict()
    resdict = result.params.valuesdict()
    print resdict['tp']

    allwav = np.logspace(np.log10(5),np.log10(1000),1024)
    src.sed.setWav(allwav)
    plt.plot(src.sed.wav, src.sed.getPL(), \
             src.sed.wav, src.sed.getBB(), \
             src.sed.wav, src.sed.getPAH(), \
             src.sed.wav, src.sed.getSED(), \
             wav, flux)
    plt.axis([6,600,0.001,2])
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
