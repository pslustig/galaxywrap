import galaxywrap as gw
import numpy as np
import pytest
from astropy.convolution import Gaussian2DKernel
from astropy.io.fits.header import Header
from astropy.wcs import WCS


def make_setup():
    setup = gw.imageproperties(1, 2, 3, 4, 5, 'ELECTRONS')
    return setup


def make_psf():
    psf = Gaussian2DKernel(x_stddev=5, x_size=50)
    psf = gw.psf(psf.array, (100, 100), 2)
    return psf


def test_setup_pixscale_scalar():
    setup = make_setup()
    setup.platescale = 2
    assert np.allclose(setup.platescale, (2, 2))


def test_mapproperties_init():
    props = make_setup()
    assert props.magzpt == 1
    assert props.exptime == 2
    assert np.allclose(props.platescale, (3, 3))
    assert props.gain == 4
    assert props.ncombine == 5
    assert props.unit == 'ELECTRONS'


def test_mapproperties_reader():
    h = Header()
    wcs = WCS()
    h['EXPTIME'] = 2
    h['GAIN'] = 4
    h['NCOMBINE'] = 5
    h['UNIT'] = 'ELECTRONS'
    h['MAGZPT'] = 1
    h.extend(wcs.to_header())

    props = gw.imageproperties.read(h)
    assert props.magzpt == 1
    assert props.exptime == 2
    assert np.allclose(props.platescale, (1, 1))
    assert props.gain == 4
    assert props.ncombine == 5
    assert props.unit == 'ELECTRONS'
