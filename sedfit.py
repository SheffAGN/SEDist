import numpy as np
import glob
from astropy import constants as const
from lmfit import minimize, Parameters
import matplotlib.pyplot as plt

class sed():

    #Speed of light in um/s:
    c = 1e6*(const.c).value
    h = (const.h).value
    kb = (const.k_B).value

    def __init__(self, \
                 norm=1., temp=25., beta=1.5, \
                 plnorm=0.5, alpha=2., turnover=50., \
                 pah=0.1):

        #Set the wavelength array; can be changed with setwav():
        self.wav = np.logspace(np.log10(5),np.log10(1000),1024)
        self.nu = self.c/self.wav

        #Set default properties; can be changed with setbb():
        self.params = {'bb_norm':norm,\
                       'bb_temp':temp,\
                       'beta':beta,\
                       'plnorm':plnorm,\
                       'alpha':alpha,\
                       'tp':turnover,\
                       'pahnorm':pah}

        #Read in the PAH features:
        self.pwav, self.pflux = np.loadtxt('M11_SEDs/Avg/avg_pah.dat', \
                                 comments='#', \
                                 unpack=True)
        self.pah = np.interp(self.wav,self.pwav,self.pflux,\
                             left=0.,right=0.)

    def setwav(self, wav):
        #(Re)Define the wavelength array (in microns):
        self.wav = wav
        self.nu = self.c/self.wav
        self.pah = np.interp(self.wav,self.pwav,self.pflux,\
                             left=0.,right=0.)

    def setbb(self, **kwargs):
        #(Re)Define the Blackbody parameters:
        #kwargs is used to maintain current values by default
        if 'bbnorm' in kwargs:
            self.params['bb_norm'] = kwargs.get('bbnorm')
        if 'temp' in kwargs:
            self.params['bb_temp'] = kwargs.get('temp')
        if 'beta' in kwargs:
            self.params['beta'] = kwargs.get('beta')

    def setpl(self, **kwargs):
        #(Re)Define the Blackbody parameters:
        if 'plnorm' in kwargs:
            self.params['plnorm'] = kwargs.get('plnorm')
        if 'alpha' in kwargs:
            self.params['alpha'] = kwargs.get('alpha')
        if 'turnover' in kwargs:
            self.params['tp'] = kwargs.get('turnover')

    def setpah(self, **kwargs):
        if 'pahnorm' in kwargs:
            self.params['pahnorm'] = kwargs.get('pahnorm')

    def setsed(self, **kwargs):
        #This is just a shorthand for setting the three SED components:
        self.setpl(**kwargs)
        self.setbb(**kwargs)
        self.setpah(**kwargs)

    def getwav(self):
        return self.wav

    def bb(self, wav, temp):

        #Use own BB generator (ripped from astropy);
        #astropy's is slow, largely due to checks.
        nu = self.c/wav
        log_boltz = self.h * nu / (self.kb * temp)
        boltzm1 = np.expm1(log_boltz)
        return (2.0 * self.h * nu ** 3 / \
                (self.c ** 2 * boltzm1))

    def getbb(self):

        #Get the modified BB spectrum:
        #Calculate normalisation:
        wavpeak = 2.9e3/self.params['bb_temp']
        nupeak = self.c/wavpeak
        norm = self.params['bb_norm']/((nupeak**self.params['beta'])*\
                       self.bb(wavpeak, self.params['bb_temp']))

        bb = self.bb(self.wav, self.params['bb_temp'])
        return norm*bb*(self.nu**self.params['beta'])

    def getpl(self):
        #Get PL spectrum:
        norm = self.params['bb_norm']*self.params['plnorm']/\
            ((self.params['tp']**self.params['alpha'])*np.exp(-1.))

        return norm*(self.wav**self.params['alpha'])*\
                             np.exp(-(self.wav/self.params['tp'])**2.)

    def getpah(self):
        return self.params['pahnorm']*self.pah

    def getsed(self):
        return self.getpl()+self.getbb()+self.getpah()

    def fitfunc(self, pars, x, y, err):

        for key, value in self.params.iteritems():
            self.params[key] = pars[key].value

        return (self.getsed() - y)/err

    def fit(self, y, err):

        pars = Parameters()
        for key, value in self.params.iteritems():
            pars.add(key, value=value)
            pars[key].set(min=0)

        pars['bb_temp'].set(min=20, max=50)
        pars['beta'].set(vary=False)
        pars['tp'].set(min=10, max=70)
        pars['alpha'].set(min=0, max=4)
        return minimize(self.fitfunc, pars, args=(self.wav, y, err))

class photsys():

    def __init__(self):
        '''
        Initialise filters and flag that they've
        not yet been interpolated to a common wavelength
        array.
        '''
        self.filter = {}
        self.interp = False

    def add_filter(self, name, wav, trans):
        self.filter[name] = {'wav':wav, 'trans':trans}
        self.interp = False

    def get_flux(self, source):

        #Shift to observed frame:
        obswav = source.restsed.getwav()*(1.+source.z)
        flux = {}

        if self.interp == False:
            for key,value in self.filter.iteritems():
                self.filter[key]['trans'] = np.interp(obswav,\
                    value['wav'], value['trans'],\
                    left=0., right=0.)
            self.interp = True

        for key,value in self.filter.iteritems():
            flux[key] = value['trans']*source.restsed.getsed()

        #plt.plot(obswav, source.restsed.getsed())
        #plt.plot(obswav, flux['red'])
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.show()

class source():

    def __init__(self, ra, dec, z):

        #Set the position and redshift of the source
        self.pos = [ra, dec]
        self.z = z

        #Give the source a rest-frame SED:
        self.restsed = sed()

#Loop through the files and fit:
a = source(0,0,1)
phot = photsys()
fwav = np.linspace(100,400,128)
ftrans = np.zeros(128)
o = (fwav > 200) & (fwav < 300)
ftrans[o] = 1.
phot.add_filter('red', fwav, ftrans)
phot.get_flux(a)

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
