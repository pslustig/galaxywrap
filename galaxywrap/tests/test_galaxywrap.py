import galaxywrap as gw
import astropy.units as u
import numpy as np
import pytest
from astropy.convolution import Gaussian2DKernel


def make_setup():
    setup = gw.setup(1, (1, 1) * u.arcsec, 'regular', 0)
    return setup


def make_psf():
    psf = Gaussian2DKernel(x_stddev=5, x_size=50)
    psf = gw.psf(psf.array, (100, 100), 2)
    return psf


def test_setup_pixscale_scalar():
    setup = make_setup()
    setup.platescale = 2 * u.arcsec
    assert u.allclose(setup.platescale, (2, 2) * u.arcsec)


def test_setup_pixscale_scalar_arcmin():
    setup = make_setup()
    setup.platescale = 2. / 60 * u.arcmin
    assert u.allclose(setup.platescale, (2, 2) * u.arcsec)


def test_displaytype_valid():
    setup = make_setup()
    validtypes = ['regular', 'curses', 'both']
    for validtype in validtypes:
        setup.displaytype = validtype
        assert setup.displaytype == validtype


def test_displaytype_invalid():
    setup = make_setup()
    with pytest.raises(AssertionError):
        setup.displaytype = 'sggseg'


def test_runoption_valid():
    setup = make_setup()
    validoptions = range(3)
    for validoption in validoptions:
        setup.runoption = validoption
        assert setup.runoption == validoption


def test_runoption_invalid():
    setup = make_setup()
    invalidoptions = [2.5, 3]
    for invalidoption in invalidoptions:
        with pytest.raises(AssertionError):
            setup.runoption = invalidoption

def test_psf_init():
    pass
