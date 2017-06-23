import numpy as np
import glob
import matplotlib.pyplot as plt
from sedfit import source
import corner

fig = plt.figure(1)
i = 1
pars = np.zeros([14,7])

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

    j = 0
    for key, value in resdict.iteritems():
        pars[i-1,j] = value
        j = j+1

    allwav = np.logspace(np.log10(5),np.log10(1000),1024)
    src.sed.setWav(allwav)

    ax = fig.add_subplot(4,4,i)
    ax.plot(src.sed.wav, src.sed.getPL(), 'b-', alpha=0.3)
    ax.plot(src.sed.wav, src.sed.getBB(), 'm-', alpha=0.3)
    ax.plot(src.sed.wav, src.sed.getPAH(), 'g-', alpha=0.3)
    ax.plot(src.sed.wav, src.sed.getSED(),'r-')
    ax.plot(wav, flux, 'k-')

    for tick in ax.yaxis.get_major_ticks():
        if i % 4 == 1:
            tick.label1On = True
        else:
            tick.label1On = False

    for tick in ax.xaxis.get_major_ticks():
        if i > 10:
            tick.label1On = True
        else:
            tick.label1On = False


    ax.axis([6,600,0.001,2])
    plt.xscale('log')
    plt.yscale('log')

    i = i+1
plt.figtext(0.5, 0.02, 'Wavelength / $\mu m$', ha='center')
plt.figtext(0.02, 0.5, 'Relative $\\nu F_\\nu$', va='center', \
            rotation='vertical')

#plt.show()
pars = np.delete(pars, (2,3), axis=1)
labels = ['PAH', 'TP', 'Temp', 'PLNorm', 'alpha']
fig = corner.corner(pars, figsize=(2,2), dpi=80,  \
                    labels=labels, plot_contours=False)
fig.savefig("/Users/James/Desktop/triangle.pdf")
