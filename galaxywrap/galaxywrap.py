from . import models
from . import utils
from astropy.nddata import NDDataArray
import numpy as np
from astropy.nddata import NDDataArray
from astropy import wcs
from astropy.io import fits
import warnings


class mapproperties(object):
    ''' setup class that contains map and run properties

    Parameters
    ----------
    fmagzeropoint : scalar
        magnitude zeropoint

    platescale : `astropy.units.quantity.Quantity`
        The angular size of one pixel. If single value the same size is
        assumed for both directions.

    '''

    # TODO: add exptime, magzpt und die beiden adu und dings faktoren

    def __init__(self, magzpt, exptime, platescale, gain, ncombine, unit):
        self.magzpt = magzpt
        self.platescale = platescale

    @property
    def platescale(self):
        return self._platescale

    @platescale.setter
    def platescale(self, platescale):

        if isinstance(platescale, wcs.WCS):
            platescale = wcs.utils.proj_plane_pixel_scales(wcs)

        if not utils.isiterable(platescale):
            platescale = [platescale, platescale]

        self._platescale = platescale

    @classmethod
    def read(cls, header, *args, **kwargs):
        if not isinstance(header, fits.header.Header):
            header = fits.getheader(header, *args, **kwargs)

        values = {}

        keys = [['EXPTIME', 'TEXPTIME'],
                ['GAIN', 'CCDGAIN'], ['NCOMBINE'], ['UNIT', 'BUNIT']]
        names = ['exptime', 'gain', 'ncombine', 'unit']

        values['platescale'] = wcs.WCS(header)
        values['magzpt'] = utils.read_zeropoint_magnitude()

        for name, key in zip(names, keys):
            values[name] = utils.read_value_or_warn(key, header)

        return cls(**values)


class psf(NDDataArray):
    # TODO: Maybe add gaussian center fit to check if psf is well centered
    def __init__(self, psf, convolutionbox, finesampling=1):
        self.data = psf
        self.convolutionbox = convolutionbox
        self.finesampling = finesampling

    @property
    def convolutionbox(self):
        return self.__convolutionbox

    @convolutionbox.setter
    def convolutionbox(self, convolutionbox):
        if not utils.isiterable(convolutionbox):
            convolutionbox = (convolutionbox, convolutionbox)
        assert isinstance(convolutionbox[0], int)
        assert isinstance(convolutionbox[1], int)
        self.__convolutionbox = convolutionbox

    @property
    def finesampling(self):
        return self.__finesampling

    @finesampling.setter
    def finesampling(self, finesampling):
        assert isinstance(finesampling, int), (
                    'finesampling factor must be an integer')
        self.__finesampling = finesampling


class image(NDDataArray):
    ''' soll maskierten array ausgeben wie nikamap, enthält settings, mask,
        uncertainty
    '''

    def __init__(self, data, *args, **kwargs):

        super(image, self).__init__(*args, **kwargs)

    def __array__(self):
        """  Overrite NDData.__array__ to force for MaskedArray output  """
        return np.ma.array(self.data, mask=self.mask)
