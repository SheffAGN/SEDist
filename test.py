import numpy as np
from sedfit import source, photset

pho = photset()

wav,trans = np.loadtxt('./sedfit/Filters/PACS/PacsFilter_red.txt', \
                       unpack=True)
pho.addFilter('P160',wav,trans)
wav,trans = np.loadtxt('./sedfit/Filters/PACS/PacsFilter_blue.txt', \
                       unpack=True)
pho.addFilter('P070',wav,trans)
wav,trans = np.loadtxt('./sedfit/Filters/PACS/PacsFilter_green.txt', \
                       unpack=True)
pho.addFilter('P100',wav,trans)

src = source(0,0,0.1)
flux = pho.getFlux(src)

ind = np.argsort(flux.keys())
print np.array(flux.keys())[ind]
print np.array(flux.values())[ind]

from sedfit import fitter
fitter.fit1d(0.,src,pho)
