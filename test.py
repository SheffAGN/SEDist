import numpy as np
from pymc3 import Model, Uniform, DensityDist
from pymc3 import NUTS, sample
from sedfit import source, photset

#Create a source:
src = source(0,0,0.1)
src.sed.setBB(temp=36.)
src.sed.setPL(alpha=3.0, turnover=40, plnorm=0.2)

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
wav,trans = np.loadtxt('./sedfit/Filters/MIPS/MipsFilter_24um.txt', \
                       unpack=True)
pho.addFilter('M024',wav,trans)
wav,trans = np.loadtxt('./sedfit/Filters/WISE/WiseFilter_W3.txt', \
                       unpack=True)
pho.addFilter('W012',wav,trans)
#wav,trans = np.loadtxt('./sedfit/Filters/WISE/WiseFilter_W2.txt', \
#                       unpack=True)
#pho.addFilter('W008',wav,trans)

obsflux = pho.getFlux(src)
sigma = 0.1*obsflux.tag.test_value
Y = obsflux+np.random.randn(5)*sigma

sedmodel = Model()

from pymc3 import find_MAP
from scipy import optimize
from pymc3 import NUTS, sample, df_summary, summary

with sedmodel:
    tp = Uniform('tp', lower=20., upper=50.)
    temp = Uniform('temp', lower=20, upper=50)
    alpha = Uniform('alpha', lower=2.7, upper=4.2)
    plnorm = Uniform('plnorm', lower=0.1, upper=0.5)
    src.sed.setBB(temp=temp)
    src.sed.setPL(alpha=alpha,turnover=tp,plnorm=plnorm)
    modflux = pho.getFlux(src)

    def logp(obs):
        return -0.5*((modflux-obs)/sigma)**2.

    Y_obs = DensityDist('Y_obs', logp, observed=Y)

    trace = sample(5000, tune=500)

    # obtain starting values via MAP
    #start = find_MAP(fmin=optimize.fmin_powell)

    # instantiate sampler
    #step = NUTS(scaling=start)

    # draw 2000 posterior samples
    #trace = sample(5000, step, start=start)

out = np.array([trace['temp'],trace['alpha'],trace['tp'], trace['plnorm']])
import corner
print df_summary(trace)
labels = ['Temp', 'alpha', 'TP', 'plnorm']
fig = corner.corner(out.T,labels=labels)
fig.savefig("/Users/James/Desktop/out.pdf")
