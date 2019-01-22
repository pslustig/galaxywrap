from galaxywrap import utils
import astropy.units as u


def test_isiterable_True():
    assert utils.isiterable((3, 4))


def test_isiterable_False():
    assert not utils.isiterable(3)


def test_change_tuple_unit():
    tpl = (5 * u.kg, 5 * u.g)
    tpl = utils.change_tuple_unit(tpl, u.kg)
    assert tpl[0].unit == u.kg
    assert tpl[1].unit == u.kg
