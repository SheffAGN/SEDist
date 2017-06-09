import numpy as np
from sedfit import source, photset
from scipy.special import erf
import matplotlib.pyplot as plt

ga = np.linspace(-2,2,1000)
si = 0.2
inerf = np.sqrt(1./2.)*(1.-(ga/si))
f = 0.5*(1.+erf(inerf))
plt.plot(ga, f, 'r-')
inerf = np.sqrt(1./2.)*(3.-(ga/si))
f = 0.5*(1.+erf(inerf))
plt.plot(ga, f, 'b-')
plt.show()
quit()


#This initialises a single source at [ra,dec]=(0,0) and z=1:
a = source(0,0,0.1)

#Whose SED properties can be adjusted thus:
a.sed.params['plnorm'] = 0.2

#This initialises a set of photometric filters:
phot = photset()

#Which you can add filters to thus:
phot.add_filter('W12', np.linspace(11,13,128), np.ones(128))

#Passing a source to the photometric set retuns the source fluxes
#through all filters in the set:
plnorms = np.random.normal(0.5, 0.1, 10000)

"""
for file in glob.glob('M11_SEDs/Indiv/*.dat'):


    wav, flux = np.loadtxt(file, \
                            comments='#', \
                            unpack=True)
    o = np.where(wav < 500.)
    wav = wav[o]
    flux = flux[o]

    flux = flux/np.max(flux)

    sed = galsed()
    sed.setwav(wav)
    sed.setpl(alpha=2.8, plnorm=0.2, turnover=40)
    sed.setbb(temp=40., beta=1.5)

    err = flux/5.
    o = np.where(wav < 35.)
    err[o] = 1.0*err[o]

    result = sed.fit(flux, err)

    result.params.pretty_print()
    print result.redchi

    allwav = np.logspace(np.log10(5),np.log10(1000),1024)
    sed.setwav(allwav)
    plt.plot(wav, flux, \
             #wav, sed.getpl(), \
             #wav, sed.getbb(), \
             #wav, sed.getpah(), \
            sed.getwav(), sed.getsed())
    plt.axis([6,600,0.001,2])
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
"""
