import galaxywrap as gw
import numpy as np
from astropy.convolution import Gaussian2DKernel
from astropy.io.fits.header import Header
from astropy.wcs import WCS
import galaxywrap.components as cp
from astropy.nddata import StdDevUncertainty


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
    assert np.allclose(props.platescale, (3600, 3600))
    assert props.gain == 4
    assert props.ncombine == 5
    assert props.unit == 'ELECTRONS'


def test_model_empty():
    m = gw.model()
    assert len(m) == 0


def make_map(mask=True, unc=True):
    map = np.random.randn(100, 100)

    m = None
    if mask:
        m = np.random.choice([True, False], size=map.shape)

    u = None
    if unc:
        u = StdDevUncertainty(np.ones(map.size))

    map = gw.image(map, mask=m, uncertainty=u, properties=make_setup())
    return map


def test_model_sersic_init():
    c0 = cp.sersic(1, 2, 3, 4, 7, 5, 6, uncertainties={'n': 7.5},
                    fixed={'n': True})
    c1 = cp.sersic(1, 2, 3, 4, 7, 5, 6, uncertainties={'n': 7.5},
                    fixed={'n': True})

    m = gw.model([c0, c1], [True, False])
    ref = (c0, c1)

    # test __len__
    assert len(m) == 2

    # test __getitem__
    assert m[0] == c0
    assert m[1] == c1

    # test __iter__
    for mm, ref in zip(m, ref):
        assert mm == ref

    # test __delitem__
    m.__delitem__(0)
    assert m[0] == c1

    # test __setitem__
    m[0] = c0
    assert m[0] == c0


def test_model_make_head():
    psf = gw.psf(np.arange(4).reshape((2, 2)), 3, 2)
    image = make_map()
    model = gw.model()
    constraints = ''
    fitarea = ((1, 40), (2, 30))
    head = model.make_head(0, image, psf, constraints, fitarea)


def test_make_galfit_fitarea():
    image = make_map()
    fitarea = None
    assert np.array_equal(((1, image.shape[1]), (1, image.shape[0])),
                         gw.model.make_galfit_fitarea(fitarea, image))

    fitarea = np.array(((1, 2), (3, 4)))
    ref = fitarea.copy() + 1
    assert np.array_equal(ref, gw.model.make_galfit_fitarea(fitarea))
