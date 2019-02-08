import galaxywrap.models as mod
import pytest
import numpy as np
import galaxywrap as gw
from astropy.convolution import Gaussian2DKernel
from astropy.nddata import StdDevUncertainty


def make_setup():
    setup = gw.imageproperties(1, 2, 3, 4, 5, 'ELECTRONS')
    return setup


def make_psf():
    psf = Gaussian2DKernel(x_stddev=5, x_size=50)
    psf = gw.psf(psf.array, 3, 2)
    return psf


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


def test_parameter_init():
    p = mod.parameter(value=3)
    assert p.value == 3

    with pytest.raises(TypeError):
        mod.parameter(uncertainty=3)

    a = mod.parameter(value=3, uncertainty=1)
    assert a.value == 3
    assert a.uncertainty == 1


def test_parameter_bounds():
    a = mod.parameter(1, bounds=(0, 3))
    assert a.bounds == (0, 3)

    with pytest.raises(AssertionError):
        mod.parameter(-1, bounds=(0, 3))

    with pytest.raises(AssertionError):
        mod.parameter(2, bounds=(0, None))


def test_parameter_rbounds():
    a = mod.parameter(1, rbounds=(0, 3))
    assert a.rbounds == (0, 3)

    a = mod.parameter(1, rbounds=1)
    assert a.rbounds == (1, 1)

    with pytest.raises(AssertionError):
        mod.parameter(2, bounds=(0, None))


def test_parameter__repr__():
    a = mod.parameter(value=3, uncertainty=1)
    assert a.__repr__() == 'parameter value: 3, uncertainty: 1'


def test_analytic_component_init():
    c = mod.analytic_component(
                        'cname', 1, 2, 3, 4, 5, 6,
                        uncertainties={'r': 4.5, 'ar': 5.5, 'pa': 6.5},
                        fixed={'r': True, 'ar': True, 'pa': True},
                        bounds={'x': (.5, 1.5)}, rbounds={})

    assert c.r.value == 4
    assert c.ar.value == 5
    assert c.pa.value == 6

    assert c.r.uncertainty == 4.5
    assert c.ar.uncertainty == 5.5
    assert c.pa.uncertainty == 6.5

    assert c.r.fixed is True
    assert c.ar.fixed is True
    assert c.pa.fixed is True


def test_sersic_init():
    c = mod.sersic(
        1, 2, 3, 4, 7, 5, 6,
        uncertainties={'n': 7.5},
        fixed={'n': True})

    assert c.n.value == 7
    assert c.n.uncertainty == 7.5
    assert c.n.fixed is True

'''
def test_sersic_to_galfit():
    c = mod.sersic(
            1, 2, 3, 4, 7, 5, 6,
            uncertainties={'n': 7.5},
            fixed={'n': True})

    assert c.to_galfit() == ' 0) sersic                    # object name\n 1)   1.0000    2.0000  0  0 # position x, y\n 2)   3.0000         1        # total magnitude\n 3)   4.0000         1        # effective radius\n 4)   7.0000         0        # sersic index\n 8)   5.0000         1        # axis ratio\n 9)   6.0000         1        # position angle\n Z) 0                         # Skip this model in output image?(yes=1, no=0)'
'''

def test_model_empty():
    m = mod.model()
    assert len(m) == 0


def test_model_sersic_init():
    c0 = mod.sersic(1, 2, 3, 4, 7, 5, 6, uncertainties={'n': 7.5},
                    fixed={'n': True})
    c1 = mod.sersic(1, 2, 3, 4, 7, 5, 6, uncertainties={'n': 7.5},
                    fixed={'n': True})

    m = mod.model([c0, c1], [True, False])
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
    model = mod.model()
    constraints = ''
    fitarea = ((1, 40), (2, 30))
    head = model.make_head(0, image, psf, constraints, fitarea)


def test_make_galfit_fitarea():
    image = make_map()
    fitarea = None
    assert np.array_equal(((1, image.shape[1]), (1, image.shape[0])),
                          mod.model.make_galfit_fitarea(fitarea, image))

    fitarea = np.array(((1, 2), (3, 4)))
    ref = fitarea.copy() + 1
    assert np.array_equal(ref, mod.model.make_galfit_fitarea(fitarea))


def test_value_or_nan_if_None():
    f = mod.parameter.value_or_nan_if_None

    assert np.isnan(f(None))
    assert f(1) == 1
    assert np.allclose(f(None, None), (np.nan, np.nan), equal_nan=True)
    assert np.allclose(f(1, None), (1, np.nan), equal_nan=True)
    assert np.allclose(f(None, 1), (np.nan, 1), equal_nan=True)
    assert np.allclose(f(1, 2), (1, 2), equal_nan=True)


test_value_or_nan_if_None()
