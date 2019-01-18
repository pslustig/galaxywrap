from astropy.modeling import Parameter
from astropy.units import Quantity


class parameter(object):
    def __init__(self, value, uncertainty=None, name='', description=''):
        self.value = value
        self.uncertainty = uncertainty
        self.name = name
        self.description = description

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


class sersic(object):
    pass


class sky(object):
    pass


class gaussian(object):
    pass


class psf(object):
    pass
