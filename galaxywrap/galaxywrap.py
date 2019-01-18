import astropy.units as u
from . import models
from . import utils

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


class psf(object):
    def __init__(self, psf, convolutionbox, finesampling=1):
        pass
