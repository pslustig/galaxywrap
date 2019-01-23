import astropy.units as u
from . import models
from . import utils
from astropy.nddata import NDDataArray
import numpy as np
from astropy.nddata import NDDataArray


class setup(object):
    ''' setup class that contains map and run properties

    Parameters
    ----------
    fmagzeropoint : scalar
        magnitude zeropoint

    platescale : `astropy.units.quantity.Quantity`
        The angular size of one pixel. If single value the same size is
        assumed for both directions.

    displaytype : `str`
        Must be either 'regular' or 'curses' or 'both'

    runoption : `int`
        Options: 0=normal run; 1,2=make model/imgblock & quit

    '''
    # TODO: Probably can remove runoption and displaytype
    # TODO: add exptime, magzpt und die beiden adu und dings faktoren
    # TODO: replace double underscores

    def __init__(self, magzeropoint, platescale, displaytype, runoption):
        self.magzeropoint = magzeropoint
        self.platescale = platescale
        self.displaytype = displaytype
        self.runoption = runoption

    @property
    def platescale(self):
        return self.__platescale

    @platescale.setter
    def platescale(self, platescale):
        assert isinstance(platescale, u.Quantity), ('platescale must have '
                                                    'angle equivalent unit')
        if not utils.isiterable(platescale):
            platescale = u.Quantity([platescale, platescale]).to(u.arcsec)

        self.__platescale = platescale

    @property
    def displaytype(self):
        return self.__displaytype

    @displaytype.setter
    def displaytype(self, displaytype):
        assert displaytype in ['regular', 'curses', 'both'], (
                    '{} no valid option for displaytype'.format(displaytype))
        self.__displaytype = displaytype

    @property
    def runoption(self):
        return self.__runoption

    @runoption.setter
    def runoption(self, runoption):
        assert runoption in range(3), (
                    '{} no valid option for runoption'.format(runoption))
        self.__runoption = runoption


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
    pass
