import numpy as np
from astropy import constants as const
from astropy.cosmology import WMAP9 as cosmo
from lmfit import minimize, Parameters
import matplotlib.pyplot as plt
import os

class source():

    def __init__(self, ra, dec, z):

        #Set the position and redshift of the source
        self.pos = [ra, dec]
        self.z = z
        self.d = (cosmo.luminosity_distance(z)).value*3.086e24
        self.d2 = 4.*np.pi*(self.d*self.d)

        #Give the source a rest-frame SED:
        self.sed = sed()

class sed():

    #Speed of light in um/s:
    c = 1e6*(const.c).value
    h = (const.h).value
    kb = (const.k_B).value

    def __init__(self, \
                 norm=1e35, temp=25., beta=1.5, \
                 plnorm=0.5, alpha=2., turnover=50., \
                 pah=0.1):

        #Set the wavelength array; can be changed with setwav():
        self.wav = np.logspace(np.log10(5),np.log10(1000),512)
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
        pwd = os.path.dirname(os.path.abspath(__file__))
        self.pwav, self.pflux = np.loadtxt(pwd+'/M11_SEDs/Avg/avg_pah.dat', \
                                 comments='#', \
                                 unpack=True)
        self.pah = np.interp(self.wav,self.pwav,self.pflux,\
                             left=0.,right=0.)

    def setWav(self, wav):
        #(Re)Define the wavelength array (in microns):
        self.wav = wav
        self.nu = self.c/self.wav
        self.pah = np.interp(self.wav,self.pwav,self.pflux,\
                             left=0.,right=0.)

    def setBB(self, **kwargs):
        #(Re)Define the Blackbody parameters:
        #kwargs is used to maintain current values by default
        if 'bbnorm' in kwargs:
            self.params['bb_norm'] = kwargs.get('bbnorm')
        if 'temp' in kwargs:
            self.params['bb_temp'] = kwargs.get('temp')
        if 'beta' in kwargs:
            self.params['beta'] = kwargs.get('beta')

    def setPL(self, **kwargs):
        #(Re)Define the Blackbody parameters:
        if 'plnorm' in kwargs:
            self.params['plnorm'] = kwargs.get('plnorm')
        if 'alpha' in kwargs:
            self.params['alpha'] = kwargs.get('alpha')
        if 'turnover' in kwargs:
            self.params['tp'] = kwargs.get('turnover')

    def setPAH(self, **kwargs):
        if 'pahnorm' in kwargs:
            self.params['pahnorm'] = kwargs.get('pahnorm')

    def setSED(self, **kwargs):
        #This is just a shorthand for setting the three SED components:
        self.setPL(**kwargs)
        self.setBB(**kwargs)
        self.setPAH(**kwargs)

    def BB(self, wav, temp):

        #Use own BB generator (ripped from astropy);
        #astropy's is slow, largely due to checks.
        nu = self.c/wav
        log_boltz = self.h * nu / (self.kb * temp)
        boltzm1 = np.expm1(log_boltz)
        return (2.0 * self.h * nu ** 3 / \
                (self.c ** 2 * boltzm1))

    def getBB(self):

        #Get the modified BB spectrum:
        #Calculate normalisation:
        wavpeak = 2.9e3/self.params['bb_temp']
        nupeak = self.c/wavpeak
        norm = self.params['bb_norm']/((nupeak**self.params['beta'])*\
                       self.BB(wavpeak, self.params['bb_temp']))

        bb = self.BB(self.wav, self.params['bb_temp'])
        return norm*bb*(self.nu**self.params['beta'])

    def getPL(self):
        #Get PL spectrum:
        norm = self.params['bb_norm']*self.params['plnorm']/\
            ((self.params['tp']**self.params['alpha'])*np.exp(-1.))

        return norm*(self.wav**self.params['alpha'])*\
                             np.exp(-(self.wav/self.params['tp'])**2.)

    def getPAH(self):
        return self.params['bb_norm']*self.params['pahnorm']*self.pah

    def getSED(self):
        return self.getPL()+self.getBB()+self.getPAH()

    def fitFunc(self, pars, x, y, err):

        for key, value in self.params.iteritems():
            self.params[key] = pars[key].value

        return (self.getSED() - y)/err

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
