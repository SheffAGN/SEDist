import numpy as np
from pymc3 import Model, Normal, HalfNormal
from pymc3 import NUTS, sample

class fitter():

    def __init__(self):
        pass

    @staticmethod
    def fit1d(data, source, photset):

        sedmodel = Model()

        with sedmodel:

            # Priors for unknown model parameters
            plnorm = Normal('plnorm', mu=0.5, sd=1.)
            print plnorm.init_value
            source.sed.setPL(plnorm=plnorm)

            sed = source.sed.getSED()
            print sed.tag.test_value.shape
