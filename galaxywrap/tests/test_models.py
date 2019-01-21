import galaxywrap.models as mod
import astropy.units as u
import pytest


def test_parameter_init_value():
    p = mod.parameter(value=3)
    assert p.value == 3
    assert p.unit is None


def test_parameter_init_novalue():
    with pytest.raises(TypeError):
        mod.parameter(uncertainty=3)


def test_parameter_init_uncertainty():
    a = mod.parameter(value=3, uncertainty=1)
    assert a.uncertainty == 1


def test_parameter_value_unit():
    p = mod.parameter(value=5*u.pix)
    assert p.value == 5 * u.pix
    assert p.unit == u.pix


def test_parameter__repr__():
    a = mod.parameter(value=3, uncertainty=1)
    assert a.__repr__() == 'parameter value: 3, uncertainty: 1'


def test_parameter_change_unc_unit():
    p = mod.parameter(5*u.kg, 500*u.g)
    assert p.uncertainty.unit == u.kg

    p.value = 3000 * u.g
    assert p.uncertainty.unit == u.g


def test_component_init():
    c = mod.component('psf', 5000, 4, 20,
                      uncertainties={'x': 1, 'y': 2, 'mag': 3},
                      fixed={'x': True, 'y': True, 'mag': True})

    assert c.x.value == 5000
    assert c.y.value == 4
    assert c.mag.value == 20

    assert c.x.uncertainty == 1
    assert c.y.uncertainty == 2
    assert c.mag.uncertainty == 3

    assert c.x.fixed is True
    assert c.y.fixed is True
    assert c.mag.fixed is True


def test_component_to_galfit():
    c = mod.component('psf', 5000, 4, 20,
                      uncertainties={'x': 1, 'y': 2, 'mag': 3},
                      fixed={'x': True, 'y': True, 'mag': True})

    assert c.to_galfit() == ''' 0) psf                       # object name\n 1) 5000.0000    4.0000  1  1 # position x, y\n 2)  20.0000         0        # total magnitude\n Z) 0                         # Skip this model in output image?(yes=1, no=0)'''


def test_analytic_component_init():
    c = mod.analytic_component(
                        'cname', 1, 2, 3, 4, 5, 6,
                        uncertainties={'r': 4.5, 'ratio': 5.5, 'pa': 6.5},
                        fixed={'r': True, 'ratio': True, 'pa': True})

    assert c.r.value == 4
    assert c.ratio.value == 5
    assert c.pa.value == 6

    assert c.r.uncertainty == 4.5
    assert c.ratio.uncertainty == 5.5
    assert c.pa.uncertainty == 6.5

    assert c.r.fixed is True
    assert c.ratio.fixed is True
    assert c.pa.fixed is True


def test_analytic_component_to_galfit():
    c = mod.analytic_component(
                        'cname', 1, 2, 3, 4, 5, 6,
                        uncertainties={'r': 4.5, 'ratio': 5.5, 'pa': 6.5},
                        fixed={'r': True, 'ratio': True, 'pa': True})
    assert c.to_galfit() == ' 0) cname                     # object name\n 1)   1.0000    2.0000  0  0 # position x, y\n 2)   3.0000         1        # total magnitude\n 3)   4.0000         0        # effective radius\n 8)   5.0000         0        # axis ratio\n 9)   6.0000         0        # position angle\n Z) 0                         # Skip this model in output image?(yes=1, no=0)'
