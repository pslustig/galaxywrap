from . import models
from . import utils
from astropy.nddata import NDDataArray, StdDevUncertainty
import numpy as np
from astropy import wcs
from astropy.io import fits
import warnings


class imageproperties(object):
    ''' setup class that contains map and run properties

    Parameters
    ----------
    fmagzeropoint : scalar
        magnitude zeropoint

    platescale : `astropy.units.quantity.Quantity`
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
            platescale = wcs.utils.proj_plane_pixel_scales(platescale)

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
    ''' soll maskierten array ausgeben wie nikamap, enthält settings, mask,
        uncertainty
    '''

    def __init__(self, data, *args, **kwargs):

        properties = kwargs.pop('properties', imageproperties())
        super(image, self).__init__(data, *args, **kwargs)
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
                uncertainty = StdDevUncertainty(hdul[uncidx].data)

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
