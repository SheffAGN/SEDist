import numpy as np
from pymc3 import Model, Normal, DensityDist
from sedfit import source, photset

#Create a source:
src = source(0,0,0.1)
src.sed.setBB(temp=40.)
src.sed.setPL(alpha=2.)

#Generate a filter set and add filters:
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

sigma = 20.
obsflux = pho.getFlux(src)
Y = obsflux.values()+np.random.randn(3)*sigma

sedmodel = Model()
with sedmodel:
    temp = Normal('temp', mu=30, sd=10)
    alpha = Normal('alpha', mu=2., sd=0.5)
    src.sed.setBB(temp=temp)
    src.sed.setPL(alpha=alpha)
    modflux = (pho.getFlux(src)).values()
    print type(modflux)
    def logp(obs):
        return -0.5*((modflux-obs)/sigma)**2.

    Y_obs = DensityDist('Y_obs', logp, observed=Y)
quit()
from pymc3 import find_MAP
from scipy import optimize
map_estimate = find_MAP(model=sedmodel, fmin=optimize.fmin_powell)
print(map_estimate)

from pymc3 import NUTS, sample
from scipy import optimize
with sedmodel:
    # draw 500 posterior samples
    trace = sample(5500)
print trace['temp']

import matplotlib.pyplot as plt
from pymc3 import traceplot, summary
traceplot(trace)
summary(trace)
