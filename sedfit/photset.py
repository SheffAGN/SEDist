import numpy as np
import theano.tensor as T

class photset():

    def __init__(self):
        '''
        Initialise filters and flag that they've
        not yet been interpolated to a common wavelength
        array.
        '''
        self.filter = {}
        self.interp = False

    def addFilter(self, name, wav, trans):
        self.filter[name] = {'wav':wav, 'trans':trans}
        self.interp = False

    def getFlux(self, source):

    #Shift to observed frame:
        obswav = source.sed.wav*(1.+source.z)
        obsnu = source.sed.c/obswav
        obssed = (source.sed.getSED()/source.d2)/1e-23

    #Interpolate filter responses onto common wavelength array:
        i = 0
        if self.interp == False:
            self.tarray = np.zeros([obswav.size, len(self.filter)])
            for key,value in self.filter.iteritems():
                o = value['trans'] > 0.
                self.tarray[:,i] = np.interp(obswav,\
                    value['wav'][o], value['trans'][o],\
                    left=0., right=0.)
                i = i+1
            self.interp = True

        #Integrate over filter transmission curves:
        obssed = obssed[:,None]

        fx = self.tarray*obssed
        gx = self.tarray

        dx = obsnu[1:] - obsnu[0:-1]
        dx = dx[:,None]

        f = dx*(fx[1:,:]+fx[0:-1,:])
        n = dx*(gx[1:,:]+gx[0:-1,:])
        norm = 1./T.sum(n, axis=0)

        flux = T.sum(f, axis=0)
        flux = norm*T.sum(f, axis=0)

        return flux
        #plt.plot(obswav, obssed)
        #plt.plot(obswav, 100.*value['trans'])
        #plt.plot(80., flux['red'], 'ro')
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.show()
