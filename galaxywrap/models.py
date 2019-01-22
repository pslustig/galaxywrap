from astropy.modeling import Parameter
from astropy.units import Quantity
from . import utils


class parameter(object):
    def __init__(self, value, uncertainty=None, name='', description='',
                 fixed=False, bounds=(None, None), rbounds=(None, None)):
        # TODO: test boundaries
        self.value = value
        self.uncertainty = uncertainty
        self.name = name
        self.description = description
        assert isinstance(fixed, bool)
        self.fixed = fixed
        self.bounds = bounds

    @property
    def uncertainty(self):
        uncertainty = self._uncertainty
        if isinstance(self.value, Quantity):
            uncertainty = uncertainty.to(self.value.unit)
        return uncertainty

    @uncertainty.setter
    def uncertainty(self, uncertainty):
        if (uncertainty is not None) and isinstance(self.value, Quantity):
                uncertainty = uncertainty.to(self.value.unit)
        self._uncertainty = uncertainty

    @property
    def unit(self):
        unit = None
        if isinstance(self.value, Quantity):
            unit = self.value.unit
        return unit

    def __repr__(self):
        return 'parameter value: {}, uncertainty: {}'.format(
                                        self.value, self.uncertainty)

    @property
    def bounds(self):
        return utils.change_tuple_unit(self._bounds, self.unit)

    @bounds.setter
    def bounds(self, bounds):
        utils.check_all_or_no_None(bounds)
        self._bounds = utils.change_tuple_unit(bounds, self.unit)

    @property
    def rbounds(self):
        return utils.change_tuple_unit(self._rbounds, self.unit)

    @rbounds.setter
    def rbounds(self, rbounds):
        utils.check_all_or_no_None(rbounds)
        self._rbounds = utils.change_tuple_unit(rbounds, self.unit)


def make_parameter(name, description, value, uncertainties, fixeddict,
                   bounddict, rbounddict):
    p = parameter(value, uncertainties.pop(name, None), name=name,
                  description=description, fixed=fixeddict.pop(name, False),
                  bounds=bounddict.pop(name, (None, None)),
                  rbounds=rbounddict.pop(name, (None, None)))
    return p


class component(object):
    # TODO: create boundary constraints
    def __init__(self, name, x, y, mag, uncertainties, fixed, bounds, rbounds):
        self.name = name

        values = [x, y, mag]
        names = ['x', 'y', 'mag']
        descriptions = ['position x', 'position y', 'total magnitude']

        for value, name, description in zip(values, names, descriptions):
            param = make_parameter(name, description, value, uncertainties,
                                   fixed, bounds, rbounds)
            self.__setattr__(name, param)

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

    def constraints_to_galfit(self, componentnumber):
        constraints = ''
        for parameter in self._parameters:
            name = utils.translate_to_constraints_names(parameter.name)
            if parameter.bounds[0] is not None:
                constraints += '\n{:10d}  {:20s}  {:8.4f} to {:8.4f}'.format(
                                      componentnumber, name, *parameter.bounds)
            if parameter.rbounds[0] is not None:
                constraints += '\n{:10d}  {:20s}  {:8.4f}  {:8.4f}'.format(
                                     componentnumber, name, *parameter.rbounds)

        return constraints


class analytic_component(component):
    def __init__(self, name, x, y, mag, r, ratio, pa, uncertainties,
                 fixed, bounds, rbounds):

        super().__init__(name, x, y, mag, uncertainties, fixed, bounds, rbounds)

        values = [r, ratio, pa]
        names = ['r', 'ratio', 'pa']
        descriptions = ['effective radius', 'axis ratio', 'position angle']

        for value, name, description in zip(values, names, descriptions):
            param = make_parameter(name, description, value, uncertainties,
                                   fixed, bounds, rbounds)
            self.__setattr__(name, param)

        self._parameters = [self.x, self.y, self.mag, self.r, None, None,
                            None, None, self.ratio, self.pa]


class sersic(analytic_component):
    def __init__(self, x, y, mag, r, n, ratio, pa, uncertainties={}, fixed={},
                 bounds={}, rbounds={}):
        super().__init__('sersic', x, y, mag, r, ratio, pa, uncertainties,
                         fixed, bounds, rbounds)

        self.n = make_parameter('n', 'sersic index', n, uncertainties, fixed,
                                bounds, rbounds)

        self._parameters = [self.x, self.y, self.mag, self.r, self.n, None,
                            None, None, self.ratio, self.pa]




class sky(object):
    pass


class gaussian(object):
    pass


class psf(object):
    pass


class model(object):
    def __init__(self, models=None, skipinimage=False):
        self.models = []
        self.skipinimage = []

        if not isinstance(models, (list, tuple)):
            models = [models]

        if not isinstance(skipinimage, (list, tuple)):
            skipinimage = [skipinimage]

        if len(skipinimage) != len(models):
            skipinimage = skipinimage * len(models)

        assert len(skipinimage) == len(models)

        for model, skip in zip(models, skipinimage):
            self.add_component(model, skip)

    def _check_data_type(self, comp, skipinimage):
        assert isinstance(comp, component)
        assert isinstance(skipinimage, bool)

    def add_component(self, comp, skipinimage=False):
        if comp is not None:
            self._check_data_type(comp, skipinimage)
            self.models.append(comp)
            self.skipinimage.append(skipinimage)

    def __len__(self):
        assert len(self.skipinimage) == len(self.models)
        return len(self.models)

    def __getitem__(self, key):
        return self.models[key]

    def __setitem__(self, key, value):
        if not isinstance(value, (list, tuple)):
            value = (value, False)

        self._check_data_type(*value)
        self.models[key] = value[0]
        self.skipinimage[key] = value[1]

    def __delitem__(self, key):
        self.models.__delitem__(key)
        self.skipinimage.__delitem__(key)

    def __iter__(self):
        return self.models.__iter__()

    def extend(self):
        pass

    def _add_skip_in_image(self, skipinimage):
        self._skipinimage.append(skipinimage)

    def fit(self, *args, **kwargs):
        '''fit model to data. input:
        map with properties, psf, fitarea

        returns
        new model with fit results and log
        '''
        return self._start_galfitrun(0, *args, **kwargs)


    def make(self, *args, **kwargs):
        return self._start_galfitrun(2, *args, **kwargs)

    def _start_galfitrun(self, mode, map, sigma, psf):
        pass
