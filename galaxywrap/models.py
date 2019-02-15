from . import utils
import numpy as np
from . import core
from astropy.table import Table, hstack, vstack


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
        # TODO:
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

        if self.bounds[0] is not None:
            constraints += '\n{:10d}  {:20s}  {:8.4f} to {:8.4f}'.format(
                              galfitcomponentnumber, name, *self.bounds)
        if self.rbounds[0] is not None:
            constraints += '\n{:10d}  {:20s}  {:8.4f}  {:8.4f}'.format(
                              galfitcomponentnumber, name, *self.rbounds)

        return value, constraints

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
        fixed = kwargs.get('fixed', {}).pop(name, False)
        bounds = kwargs.get('bounds', {}).pop(name, (None, None))
        rbounds = kwargs.get('rbounds', {}).pop(name, (None, None))

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

        self._parameters = [self.x, self.y, self.mag, self.r, self.n, None,
                            None, None, self.ar, self.pa]


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
    def __init__(self, x, y, mag, fwhm, n, ar, pa, **kwargs):
        super(sersic, self).__init__('gaussian', x, y, mag, fwhm, ar, pa,
                                     **kwargs)

        self.__delattr__('r')
        self.n = self.make_parameter('fwhm', 'FWHM', fwhm, **kwargs)

        self._parameters = [self.x, self.y, self.mag, self.fwhm, self.n, None,
                            None, None, self.ar, self.pa]


class psf(object):
    pass


class model(object):
    def __init__(self, components=None, skipinimage=False):
        self.components = []
        self.skipinimage = []

        if not isinstance(components, (list, tuple)):
            components = [components]

        if not isinstance(skipinimage, (list, tuple)):
            skipinimage = [skipinimage]

        if len(skipinimage) != len(components):
            skipinimage = skipinimage * len(components)

        assert len(skipinimage) == len(components)

        for component, skip in zip(components, skipinimage):
            self.add_component(component, skip)

    def _check_data_type(self, comp, skipinimage):
        assert isinstance(comp, component)
        assert isinstance(skipinimage, bool)

    def add_component(self, comp, skipinimage=False):
        if comp is not None:
            self._check_data_type(comp, skipinimage)
            self.components.append(comp)
            self.skipinimage.append(skipinimage)

    def __repr__(self):
        return 'galfit model containing {} component(s)'.format(len(self))

    def __len__(self):
        assert len(self.skipinimage) == len(self.components)
        return len(self.components)

    def __getitem__(self, key):
        return self.components[key]

    def __setitem__(self, key, value):
        if not isinstance(value, (list, tuple)):
            value = (value, False)

        self._check_data_type(*value)
        self.components[key] = value[0]
        self.skipinimage[key] = value[1]

    def __delitem__(self, key):
        self.components.__delitem__(key)
        self.skipinimage.__delitem__(key)

    def __iter__(self):
        return self.components.__iter__()

    def extend(self, *models):
        for model in models:
            for component, skipinimage in zip(model.components,
                                              model.skipinimage):
                self.append(component, skipinimage)

    def _add_skip_in_image(self, skipinimage):
        self._skipinimage.append(skipinimage)

    def fit(self, image, psf, gconstraints, fitarea=None, **kwargs):
        '''fit model to data. input:
        map with properties, psf, fitarea

        returns
        new model with fit results and log
        '''
        return self._start_galfitrun(0, image, psf, gconstraints, fitarea, **kwargs)

    def make(self, image, psf, **kwargs):
        return self._start_galfitrun(2, image, psf)

    def _make_global_constraints(self, gconstraints):
        return ''

    @staticmethod
    def make_galfit_fitarea(fitarea=None, image=None):
        if fitarea is None:
            ylen, xlen = image.shape
            fitarea = np.array(((0, xlen-1), (0, ylen-1)))

        fitarea = np.array(fitarea) + 1
        return fitarea

    def _start_galfitrun(self, mode, image, psf, gconstraints, fitarea,
                         **kwargs):
        entries = ''
        constraints = self._make_global_constraints(gconstraints)

        for i, (component, skip) in enumerate(zip(self, self.skipinimage)):
            entry, constraint = component._to_galfit(i)
            entries += '\n\n' + entry
            constraints += constraint

        head = self.make_head(mode, image, psf, constraints, fitarea)
        feedme = head + '\n' + entries

        return core.fit(feedme, image, psf, constraints, **kwargs)

    @staticmethod
    def make_head(runoption, image, psf, constraints, fitarea):

        convbox = (0, 0)
        if psf is not None:
            convbox = psf.convolutionbox

        fitarea = model.make_galfit_fitarea(fitarea, image)
        unc = image.uncertainty
        head = '# IMAGE and GALFIT CONTROL PARAMETERS'
        head += '\nA) inimg.fits'
        head += '\nB) imgblock.fits'
        head += '\nC) {}'.format('sigma.fits' if unc is not None else 'none')
        head += '\nD) {}'.format('psf.fits' if psf is not None else 'none')
        head += '\nE) {}'.format(1 if psf is None else psf.finesampling)
        head += '\nF) {}'.format('none' if image.mask is None else 'mask.fits')
        head += '\nG) {}'.format('none' if constraints == '' else 'mask.fits')
        head += '\nH) {}   {}   {}   {}'.format(*fitarea[0], *fitarea[1])
        head += '\nI) {}  {}'.format(*convbox)
        head += '\nJ) {}'.format(image.properties.magzpt)
        head += '\nK) {}  {}'.format(*image.properties.platescale)
        head += '\nO) {}'.format('regular')
        head += '\nP) {}'.format(runoption)
        return head

    def _to_table(self, **kwargs):
        # maybe add index
        t = Table()
        for i, comp in enumerate(self):
            componenttable = comp._to_table(**kwargs)
            componenttable['skipinimage'] = self.skipinimage[i]
            t = vstack([t, componenttable])

        return t

    def write(self, *args, **kwargs):
        # ignore warnings...
        self._to_table().write(*args, **kwargs)
