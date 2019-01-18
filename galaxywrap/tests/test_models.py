import galaxywrap.models as mod
import astropy.units as u
import pytest


def test_parameter_init_value():
    p = mod.parameter(value=3)
    assert p.value == 3


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
