import numpy.random as rand
import corner
import matplotlib.pyplot as plt
import numpy as np
from sedfit import source, photset

pho = photset()
wav,trans = np.loadtxt('./sedfit/Filters/IRAS/IRAS_12um_Filter_Resp.dat', \
                       unpack=True, usecols=(0,3))
pho.addFilter('I12',wav,trans)

wav,trans = np.loadtxt('./sedfit/Filters/IRAS/IRAS_25um_Filter_Resp.dat', \
                       unpack=True, usecols=(0,3))
pho.addFilter('I25',wav,trans)

wav,trans = np.loadtxt('./sedfit/Filters/IRAS/IRAS_60um_Filter_Resp.dat', \
                       unpack=True, usecols=(0,3))
pho.addFilter('I60',wav,trans)

wav,trans = np.loadtxt('./sedfit/Filters/IRAS/IRAS_100um_Filter_Resp.dat', \
                       unpack=True, usecols=(0,3))
pho.addFilter('I100',wav,trans)

src = source(0,0,0.01)

nsed = 1000
tp = rand.uniform(20., 50., nsed)
temp = 53.33-(2./3.)*tp+rand.normal(0., 1.5, nsed)
plnorm = rand.uniform(0.1,0.5,nsed)
alpha = rand.uniform(2.7,4.2,nsed)
pah = rand.uniform(0.02, 0.1,nsed)

src.sed.params['bb_norm'] = 1e35
x = np.zeros(nsed)
y = np.zeros(nsed)
for i in range(nsed):
    src.sed.params['tp'] = tp[i]
    src.sed.params['bb_temp'] = temp[i]
    src.sed.params['plnorm'] = plnorm[i]
    src.sed.params['alpha'] = alpha[i]
    src.sed.params['pahnorm'] = pah[i]
    flux = pho.getFlux(src)
    x[i] = flux.tag.test_value[2]/flux.tag.test_value[3]
    y[i] = flux.tag.test_value[0]/flux.tag.test_value[1]
    print i

plt.plot(x, y, 'ro')
plt.xscale('log')
plt.yscale('log')
plt.show()
