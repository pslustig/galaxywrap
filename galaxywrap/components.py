from . import utils
import numpy as np
from astropy.table import Table, hstack


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
        self.rbounds = rbounds

    def __repr__(self):
        return 'parameter value: {}, uncertainty: {}'.format(
                                        self.value, self.uncertainty)

    @property
    def bounds(self):
        return self._bounds

    @bounds.setter
    def bounds(self, bounds):
        utils.check_all_or_no_None(bounds)
        if bounds[0] is not None:
            assert self.value >= bounds[0]
            assert self.value <= bounds[1]
        self._bounds = bounds

    @property
    def rbounds(self):
        return self._rbounds

    @rbounds.setter
    def rbounds(self, rbounds):
        if not utils.isiterable(rbounds):
            rbounds = (rbounds, rbounds)
        utils.check_all_or_no_None(rbounds)
        self._rbounds = rbounds

    def _to_galfit(self, galfitcomponentnumber, parameternumber):

        value = '\n{:2d}) {:8.4f}         {:d}        # {}'.format(
                parameternumber, self.value, not self.fixed, self.description)

        constraints = ''
        name = utils.translate_to_constraints_names(self.name)

        constraints += self.make_range_constraints(
                        galfitcomponentnumber, name, self.bounds, False)
        constraints += self.make_range_constraints(
                        galfitcomponentnumber, name, self.rbounds, True)

        return value, constraints

    @staticmethod
    def make_range_constraints(galfitcomponentnumber, name, bounds,
                               isrelative):
        constraints = ''
        if bounds[0] is not None:
            constraints = '{}  {}  '.format(galfitcomponentnumber, name)
            constraints += '{0:8.4f} {2} {1:8.4f}'.format(
                                    *bounds, 'to' if not isrelative else '')
        return constraints

    @classmethod
    def value_or_nan_if_None(self, *values):
        values = np.array(values, dtype=float)
        for i, value in enumerate(values):
            if value is None:
                values[i] = np.nan

        if len(values) == 1:
            values = values[0]

        return values

    def _to_table(self, constraints=True):
        '''
        Save parameter values in table. Used to save models to table.
        '''
        unc = self.value_or_nan_if_None(self.uncertainty)
        name = self.name

        t = Table()
        t[name] = [self.value]
        t['{}_unc'.format(name)] = [unc]

        if constraints:
            leftrbound, rightrbound = self.value_or_nan_if_None(self.rbounds)
            leftbound, rightrbound = self.value_or_nan_if_None(self.bounds)

            t['{}_fixed'.format(name)] = [self.fixed]
            t['{}_left_rbound'.format(name)] = [leftrbound]
            t['{}_right_rbound'.format(name)] = [rightrbound]
            t['{}_left_bound'.format(name)] = [leftbound]
            t['{}_right_bound'.format(name)] = [rightrbound]

        return t


class global_constraint(Table):
    def __init__(self, *args, **kwargs):
        super(Table, self).__init__(*args, **kwargs)

    def add_constraint(self, component, parameter, ctype):
        pass

    def _to_galfit(self):
        pass


class component(object):
    def __init__(self, name, values, names, descriptions, **kwargs):

        self.name = name

        for value, name, description in zip(values, names, descriptions):
            param = self.make_parameter(name, description, value, **kwargs)
            self.__setattr__(name, param)

        self._parameters = []

    def __repr__(self):
        return self._to_galfit(0)[0]

    @staticmethod
    def make_parameter(name, description, value, **kwargs):

        uncertainties = kwargs.get('uncertainties', {})
        fixed = kwargs.get('fixed', {}).get(name, False)
        bounds = kwargs.get('bounds', {}).get(name, (None, None))
        rbounds = kwargs.get('rbounds', {}).get(name, (None, None))

        p = parameter(value, uncertainties.pop(name, None), name=name,
                      description=description, fixed=fixed, bounds=bounds,
                      rbounds=rbounds)
        return p

    def _to_galfit(self, componentnumber, skipinimage=False):
        componentnumber += 1    # galfit starts counting with 1
        constraints = ''

        body = '# Component number: {}'.format(componentnumber)
        body += '\n 0) {:25s} # object name'.format(self.name)
        body += '\n 1) {:8.4f}  {:8.4f}  {:d}  {:d} # position x, y'.format(
            self.x.value+1, self.y.value+1, not self.x.fixed, not self.y.fixed)

        for i, parameter in enumerate(self._parameters):
            if parameter is not None:
                pvalue, pconstr = parameter._to_galfit(componentnumber, i+1)
                constraints += pconstr
                if i > 1:
                    body += pvalue

        body += '\n Z) {:d}                         # {}'.format(
                skipinimage, 'Skip this model in output image?(yes=1, no=0)')

        return body, constraints

    def _to_table(self, **kwargs):
        '''
        Create a single row table containing this components
        '''
        t = Table()
        for parameter in self._parameters:
            if parameter is not None:
                t = hstack([t, parameter._to_table(**kwargs)])
        return t


class analytic_component(component):
    def __init__(self, name, x, y, mag, r, ar, pa, **kwargs):

        values = [x, y, mag, r, ar, pa]
        names = ['x', 'y', 'mag', 'r', 'ar', 'pa']
        descriptions = ['position x', 'position y', 'absolute magnitude',
                        'effective radius', 'axis ratio', 'position angle']

        super(analytic_component, self).__init__(
                            name, values, names, descriptions, **kwargs)

        self._parameters = [self.x, self.y, self.mag, self.r, None, None,
                            None, None, self.ar, self.pa]


class sersic(analytic_component):
    def __init__(self, x, y, mag, r, n, ar, pa, **kwargs):
        super(sersic, self).__init__('sersic', x, y, mag, r, ar, pa,
                                     **kwargs)

        self.n = self.make_parameter('n', 'sersic index', n, **kwargs)

        self._parameters[4] = self.n


class sky(component):
    def __init__(self, bkg, dbkg_dx, dbkg_dy, **kwargs):

        values = bkg, dbkg_dx, dbkg_dy
        names = 'bkg', 'dbkg_dx', 'dbkg_dy'
        descriptions = ('bkg value at center of fitting region [ADUs]',
                        'dbkg / dx (bkg gradient in x)',
                        'dbkg / dy (bkg gradient in y)')

        super(sky, self).__init__('sky', values, names, descriptions, **kwargs)
        self._parameters = [self.bkg, self.dbkg_dx, self.dbkg_dy]

    def _to_galfit(self, componentnumber, skipinimage=False):
        componentnumber += 1    # galfit starts counting with 1
        constraints = ''

        body = '# Component number: {}'.format(componentnumber)
        body += '\n 0) {:25s} # object name'.format(self.name)

        for i, parameter in enumerate(self._parameters):
            pvalue, pconstr = parameter._to_galfit(componentnumber, i+1)
            constraints += pconstr
            body += pvalue

        body += '\n Z) {:d}                         # {}'.format(
                skipinimage, 'Skip this model in output image?(yes=1, no=0)')

        return body, constraints


class gaussian(analytic_component):
    def __init__(self, x, y, mag, fwhm, ar, pa, **kwargs):
        super(sersic, self).__init__('gaussian', x, y, mag, fwhm, ar, pa,
                                     **kwargs)

        self.fwhm = self.make_parameter('fwhm', 'FWHM', fwhm, **kwargs)

        self._parameters[3] = self.fwhm
        self.__delattr__('r')


class psf(component):
    def __init__(self, x, y, mag, **kwargs):

        values = x, y, mag
        names = 'x', 'y', 'mag'
        descriptions = ('position x', 'position y', 'absolute magnitude')

        super(psf, self).__init__('psf', values, names, descriptions, **kwargs)
        self._parameters = [self.x, self.y, self.mag]
