import ROOT
import rhalphalib as rl
import numpy as np

class TransferFactorSample(rl.ParametericSample):
    def __init__(self, samplename, sampletype, transferfactor, dependentsample, nominal_values=None, stat_unc=None, observable=None, min_val=None, epsilon=1e-5, effect_threshold=0.01, channel_name=None):
        """
        Create a sample that depends on another Sample by some transfer factor.
        The transfor factor can be a constant, an array of parameters of same length
        as the dependent sample binning, or a matrix of parameters where the second
        dimension matches the sample binning, i.e. expectation = tf @ dependent_expectation.
        The latter requires an additional observable argument to specify the definition of the first dimension.
        In all cases, please use numpy object arrays of Parameter types.
        Passing in a ``min_val`` means param values will be clipped at the min_val.
        """
        if not isinstance(transferfactor, np.ndarray):
            raise ValueError("Transfer factor is not a numpy array")
        if not isinstance(dependentsample, rl.Sample):
            raise ValueError("Dependent sample does not inherit from Sample")
        if len(transferfactor.shape) == 2:
            if observable is None:
                raise ValueError("Transfer factor is 2D array, please provide an observable")
            params = np.dot(transferfactor, dependentsample.getExpectation())
            if min_val is not None:
                for idx, p in np.ndenumerate(params):
                    params[idx] = p.max(min_val)
        elif len(transferfactor.shape) <= 1:
            observable = dependentsample.observable
            if stat_unc is not None:
                name = samplename if channel_name is None else channel_name
                MCStatTemplate = (np.ones_like(stat_unc), observable._binning, observable._name)
                MCStat = rl.TemplateSample(name+'_mcstat', rl.Sample.BACKGROUND, MCStatTemplate)
                for i in range(MCStat.observable.nbins):
                    effect_up = np.ones_like(MCStat._nominal)
                    effect_down = np.ones_like(MCStat._nominal)
                    if stat_unc[i] < effect_threshold:
                        continue
                    effect_up[i] = 1.0 + min(1.0, stat_unc[i])
                    effect_down[i] = max(epsilon, 1.0 - min(1.0, stat_unc[i]))
                    print(samplename, i, nominal_values[i], effect_up[i], effect_down[i])
                    param = rl.NuisanceParameter(name + '_mcstat_bin%i' % i, combinePrior='shape')
                    MCStat.setParamEffect(param, effect_up, effect_down)
                params = transferfactor * MCStat.getExpectation() * dependentsample.getExpectation()
            else:
                params = transferfactor * dependentsample.getExpectation()
            if min_val is not None:
                for i, p in enumerate(params):
                    params[i] = p.max(min_val)
        else:
            raise ValueError("Transfer factor has invalid dimension")
        super(TransferFactorSample, self).__init__(samplename, sampletype, observable, params)
        self._transferfactor = transferfactor
        self._dependentsample = dependentsample
        self._stat_unc = stat_unc
        self._nominal_values = nominal_values

    @property
    def transferfactor(self):
        return self._transferfactor

    @property
    def dependentsample(self):
        return self._dependentsample

    @property
    def stat_unc(self):
        return self._stat_unc

    @property
    def nominal_values(self):
        return self._nominal_values
