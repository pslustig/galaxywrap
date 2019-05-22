from . import utils
from . import components as cp
from astropy.nddata import NDDataArray, StdDevUncertainty
import numpy as np
from astropy import wcs
from astropy.io import fits
import warnings
from astropy.table import Table, vstack
from . import core


class imageproperties(object):
    ''' setup class that contains map and run properties

    Parameters
    ----------
    fmagzeropoint : scalar
        magnitude zeropoint

    platescale : ``
        The angular size of one pixel. If single value the same size is
        assumed for both directions.

    '''

    def __init__(self, magzpt=None, exptime=None, platescale=None, gain=None,
                 ncombine=None, unit=None):
        self.magzpt = magzpt
        self.exptime = exptime
        self.platescale = platescale
        self.gain = gain
        self.ncombine = ncombine
        self.unit = unit

    @property
    def platescale(self):
        return self._platescale

    @platescale.setter
    def platescale(self, platescale):

        if isinstance(platescale, wcs.WCS):
            platescale = wcs.utils.proj_plane_pixel_scales(platescale) * 3600

        if not utils.isiterable(platescale):
            platescale = [platescale, platescale]

        self._platescale = platescale

    def __repr__(self):
        items = self.__dict__.items()
        out = '\n'.join('{}: {}'.format(key.replace('_', ''), value) for key,
                        value in items)
        return out

    @classmethod
    def read(cls, header, *args, **kwargs):
        if not isinstance(header, fits.header.Header):
            header = fits.getheader(header, *args, **kwargs)

        values = {}

        keys = [['EXPTIME', 'TEXPTIME'],
                ['GAIN', 'CCDGAIN'], ['NCOMBINE'], ['UNIT', 'BUNIT']]
        names = ['exptime', 'gain', 'ncombine', 'unit']

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', wcs.FITSFixedWarning)
            values['platescale'] = wcs.WCS(header)

        values['magzpt'] = utils.read_zeropoint_magnitude(header)

        for name, key in zip(names, keys):
            values[name] = utils.read_value_or_warn(key, header)

        return cls(**values)

    def to_header(self):
        header = fits.header.Header()
        header['EXPTIME'] = self.exptime
        header['UNIT'] = self.unit
        header['NCOMBINE'] = self.ncombine

        return header


class image(NDDataArray):
    ''' soll maskierten array ausgeben, enthält settings, mask,
        uncertainty
    '''

    def __init__(self, data, *args, **kwargs):

        properties = kwargs.pop('properties', imageproperties())

        unc = kwargs.pop('uncertainty', None)
        if unc is not None:
            unc = StdDevUncertainty(unc)

        super(image, self).__init__(data, uncertainty=unc, *args, **kwargs)
        self.properties = properties

    def __array__(self):
        """  Overrite NDData.__array__ to force for MaskedArray output  """
        return np.ma.array(self.data, mask=self.mask)

    @classmethod
    def read(cls, filename, dataidx=0, headeridx=None, maskidx=None,
             uncidx=None):

        mask = None
        uncertainty = None

        if headeridx is None:
            headeridx = dataidx

        with fits.open(filename) as hdul:
            data = hdul[dataidx].data
            imgprops = imageproperties.read(
                    fits.header.Header(hdul[headeridx].header))
            if maskidx is not None:
                mask = hdul[maskidx].data.astype(bool)

            if uncidx is not None:
                uncertainty = hdul[uncidx].data

        return cls(data, properties=imgprops, mask=mask,
                   uncertainty=uncertainty)


class psf(NDDataArray):
    def __init__(self, data, convolutionbox, finesampling=1):
        super(psf, self).__init__(data)
        self.convolutionbox = convolutionbox
        self.finesampling = finesampling

    @property
    def convolutionbox(self):
        return self._convolutionbox

    @convolutionbox.setter
    def convolutionbox(self, convolutionbox):
        if not utils.isiterable(convolutionbox):
            convolutionbox = (convolutionbox, convolutionbox)
        assert isinstance(convolutionbox[0], int)
        assert isinstance(convolutionbox[1], int)
        self._convolutionbox = convolutionbox

    @property
    def finesampling(self):
        return self._finesampling

    @finesampling.setter
    def finesampling(self, finesampling):
        assert isinstance(finesampling, int), (
                    'finesampling factor must be an integer')
        self._finesampling = finesampling


class model(object):
    def __init__(self, components=None, skipinimage=False):
        self.components = []
        self.skipinimage = []

        if not isinstance(components, (list, tuple)):
            components = [components]

        if isinstance(skipinimage, bool):
            skipinimage = [skipinimage]

        if len(skipinimage) != len(components):
            skipinimage = skipinimage * len(components)

        assert len(skipinimage) == len(components)

        for component, skip in zip(components, skipinimage):
            self.add_component(component, skip)

    def _check_data_type(self, comp, skipinimage):
        assert isinstance(comp, cp.component)
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

    def fit(self, image, psf, fitarea=None, gconstraints=None, **kwargs):
        '''fit model to data. input:
        map with properties, psf, fitarea

        returns
        new model with fit results and log
        '''
        return self._start_galfitrun(
                    0, image, psf, fitarea, gconstraints, **kwargs)

    def make(self, image, psf, **kwargs):
        return self._start_galfitrun(
                    2, image, psf, fitarea=None, gconstraints=None, **kwargs)

    def _make_global_constraints(self, gconstraints):
        return ''

    @staticmethod
    def make_galfit_fitarea(fitarea=None, image=None):
        if fitarea is None:
            ylen, xlen = image.shape
            fitarea = np.array(((0, xlen-1), (0, ylen-1)))

        fitarea = np.array(fitarea) + 1
        return fitarea

    def _start_galfitrun(self, mode, image, psf, fitarea, gconstraints,
                         **kwargs):
        entries = ''     # component entries
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
        head += '\nG) {}'.format(
                        'none' if constraints == '' else 'constraints.txt ')
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
