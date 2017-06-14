import numpy as np

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
        flux = {}

        #Interpolate filter responses onto common wavelength array:
        if self.interp == False:
            for key,value in self.filter.iteritems():
                o = value['trans'] > 0.
                self.filter[key]['trans'] = np.interp(obswav,\
                    value['wav'][o], value['trans'][o],\
                    left=0., right=0.)
            self.interp = True

        #Integrate over filter transmission curve:
        for key,value in self.filter.iteritems():
            fx = value['trans']*obssed
            gx = value['trans']
            dx = obsnu[1:] - obsnu[0:-1]
            norm = 1./np.sum(dx*(gx[1:]+gx[0:-1]))
            flux[key] = (norm*dx*(fx[1:]+fx[0:-1])).sum()

        return flux
        #plt.plot(obswav, obssed)
        #plt.plot(obswav, 100.*value['trans'])
        #plt.plot(80., flux['red'], 'ro')
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.show()
