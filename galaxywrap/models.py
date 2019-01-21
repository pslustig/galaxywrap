from astropy.modeling import Parameter
from astropy.units import Quantity


class parameter(object):
    def __init__(self, value, uncertainty=None, name='', description='',
                 fixed=False):
        self.value = value
        self.uncertainty = uncertainty
        self.name = name
        self.description = description
        assert isinstance(fixed, bool)
        self.fixed = fixed

    @property
    def uncertainty(self):
        uncertainty = self.__uncertainty
        if isinstance(self.value, Quantity):
            uncertainty = uncertainty.to(self.value.unit)
        return uncertainty

    @uncertainty.setter
    def uncertainty(self, uncertainty):
        if (uncertainty is not None) and isinstance(self.value, Quantity):
                uncertainty = uncertainty.to(self.value.unit)
        self.__uncertainty = uncertainty

    @property
    def unit(self):
        unit = None
        if isinstance(self.value, Quantity):
            unit = self.value.unit
        return unit

    def __repr__(self):
        return 'parameter value: {}, uncertainty: {}'.format(
                                        self.value, self.uncertainty)


class component(object):
    def __init__(self, name, x, y, mag, uncertainties, fixed):
        self.name = name
        self.x = parameter(x, uncertainties.pop('x', None), name='x',
                           description='position x',
                           fixed=fixed.pop('x', False))

        self.y = parameter(y, uncertainties.pop('y', None), name='y',
                           description='position y',
                           fixed=fixed.pop('y', False))

        self.mag = parameter(mag, uncertainties.pop('mag', None), name='mag',
                             description='total magnitude',
                             fixed=fixed.pop('mag', False))

        self._parameters = [self.x, self.y, self.mag]

    def to_galfit(self, skipinimage=False):
        out = ' 0) {:25s} # object name'.format(self.name)
        out += '\n 1) {:8.4f}  {:8.4f}  {:d}  {:d} # position x, y'.format(
                        self.x.value, self.y.value, self.x.fixed, self.y.fixed)

        for i, parameter in enumerate(self._parameters[2:]):
            if parameter is not None:
                out += '\n{:2d}) {:8.4f}         {:d}        # {}'.format(
                            i+2, parameter.value, not parameter.fixed,
                            parameter.description)
        out += '\n Z) {:d}                         # {}'.format(
                skipinimage, 'Skip this model in output image?(yes=1, no=0)')
        return out


class analytic_component(component):
    def __init__(self, name, x, y, mag, r, ratio, pa, uncertainties,
                 fixed):

        super().__init__(name, x, y, mag, uncertainties, fixed)

        self.r = parameter(r, uncertainties.pop('r'), name='r',
                           description='effective radius',
                           fixed=fixed.pop('r', False))

        self.ratio = parameter(ratio, uncertainties.pop('ratio', None),
                               name='ratio', description='axis ratio',
                               fixed=fixed.pop('ratio', False))

        self.pa = parameter(pa, uncertainties.pop('pa', None), name='pa',
                            description='position angle',
                            fixed=fixed.pop('pa', False))

        self._parameters = [self.x, self.y, self.mag, self.r, None, None,
                            None, None, self.ratio, self.pa]


class sersic(analytic_component):
    def __init__(self, x, y, mag, r, n, ratio, pa, uncertainties={}, fixed={}):
        super().__init__('sersic', x, y, mag, r, ratio, pa, uncertainties,
                         fixed)

        self.n = parameter(n, uncertainties.pop('n', None), name='n',
                           description='sersic index',
                           fixed=fixed.pop('n', False))

        self._parameters = [self.x, self.y, self.mag, self.r, self.n, None,
                            None, None, self.ratio, self.pa]




class sky(object):
    pass


class gaussian(object):
    pass


class psf(object):
    pass
