import galaxywrap.components as cp
import pytest
import numpy as np
import galaxywrap as gw
from astropy.convolution import Gaussian2DKernel
from astropy.nddata import StdDevUncertainty


def test_make_range_constraints():

    cstr = cp.parameter._make_range_constraints(1, "x", (3, 4), False)
    assert cstr.split('  ') == ['1',  'x',  '3.0000', 'to', '4.0000\n']

    cstr = cp.parameter._make_range_constraints(1, "x", (3, 4), True)
    assert cstr.split('  ') == ['1',  'x',  '3.0000', '4.0000\n']

    cstr = cp.parameter._make_range_constraints(1, "x", (None, None), True)
    assert cstr == ''


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
    p = cp.parameter(value=3)
    assert p.value == 3

    with pytest.raises(TypeError):
        cp.parameter(uncertainty=3)

    a = cp.parameter(value=3, uncertainty=1)
    assert a.value == 3
    assert a.uncertainty == 1


def test_parameter_bounds():
    a = cp.parameter(1, bounds=(0, 3))
    assert a.bounds == (0, 3)

    with pytest.raises(AssertionError):
        cp.parameter(-1, bounds=(0, 3))

    with pytest.raises(AssertionError):
        cp.parameter(2, bounds=(0, None))


def test_parameter_rbounds():
    a = cp.parameter(1, rbounds=(0, 3))
    assert a.rbounds == (0, 3)

    a = cp.parameter(1, rbounds=1)
    assert a.rbounds == (-1, 1)

    with pytest.raises(AssertionError):
        cp.parameter(2, bounds=(0, None))


def test_parameter__repr__():
    a = cp.parameter(value=3, uncertainty=1)
    assert a.__repr__() == 'parameter value: 3, uncertainty: 1'


def test_analytic_component_init():
    c = cp.analytic_component(
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
    c = cp.sersic(
        1, 2, 3, 4, 7, 5, 6,
        uncertainties={'n': 7.5},
        fixed={'n': True})

    assert c.n.value == 7
    assert c.n.uncertainty == 7.5
    assert c.n.fixed is True

'''
def test_sersic_to_galfit():
    c = cp.sersic(
            1, 2, 3, 4, 7, 5, 6,
            uncertainties={'n': 7.5},
            fixed={'n': True})

    assert c.to_galfit() == ' 0) sersic                    # object name\n 1)   1.0000    2.0000  0  0 # position x, y\n 2)   3.0000         1        # total magnitude\n 3)   4.0000         1        # effective radius\n 4)   7.0000         0        # sersic index\n 8)   5.0000         1        # axis ratio\n 9)   6.0000         1        # position angle\n Z) 0                         # Skip this model in output image?(yes=1, no=0)'
'''



def test_value_or_nan_if_None():
    f = cp.parameter.value_or_nan_if_None

    assert np.isnan(f(None))
    assert f(1) == 1
    assert np.allclose(f(None, None), (np.nan, np.nan), equal_nan=True)
    assert np.allclose(f(1, None), (1, np.nan), equal_nan=True)
    assert np.allclose(f(None, 1), (np.nan, 1), equal_nan=True)
    assert np.allclose(f(1, 2), (1, 2), equal_nan=True)
