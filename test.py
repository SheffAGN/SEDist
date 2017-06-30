import numpy as np
from pymc3 import Model, Uniform, DensityDist, Normal
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
from pymc3 import NUTS, sample, df_summary, summary, Metropolis

with sedmodel:
    tp = 35.+15.*Normal('tp', 0., 0.5)
    #temp = 35.+15.*Uniform('temp', lower=-1, upper=1)
    #alpha = 3.45+0.75*Uniform('alpha', lower=-1, upper=1)
    plnorm = 0.3+0.2*Normal('plnorm', 0., 0.5)

    #src.sed.setBB(temp=temp)
    src.sed.setPL(turnover=tp,plnorm=plnorm)
    modflux = pho.getFlux(src)

    def logp(obs):
        return -0.5*((modflux-obs)/sigma)**2.

    Y_obs = DensityDist('Y_obs', logp, observed=Y)

    trace = sample(1000, n_init=50000)

    # obtain starting values via MAP
    #start = find_MAP(fmin=optimize.fmin_powell)

    # instantiate sampler
    #step = NUTS(scaling=start)

    # draw 2000 posterior samples
    #trace = sample(5000, step, start=start)

out = np.array([35.+15.*trace['tp'], 0.3+0.2*trace['plnorm']])
import corner
print df_summary(trace)
labels = ['TP', 'plnorm']
fig = corner.corner(out.T,labels=labels, plot_density=False, plot_contours=False)
fig.savefig("out.pdf")
