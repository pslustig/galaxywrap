from galaxywrap import utils
import astropy.units as u
import pytest
from astropy.io.fits.header import Header
import numpy as np
import warnings
from astropy.io.fits.verify import VerifyWarning


def test_isiterable_True():
    assert utils.isiterable((3, 4))


def test_isiterable_False():
    assert not utils.isiterable(3)


def test_change_tuple_unit():
    tpl = (5 * u.kg, 5 * u.g)
    tpl = utils.change_tuple_unit(tpl, u.kg)
    assert tpl[0].unit == u.kg
    assert tpl[1].unit == u.kg


def test_translate_to_constraints_names():
    assert utils.translate_to_constraints_names('r') == 're'
    assert utils.translate_to_constraints_names('a') == 'a'


def test_all_or_no_None():
    utils.check_all_or_no_None((None, None))
    utils.check_all_or_no_None((1, 1))

    with pytest.raises(AssertionError):
        utils.check_all_or_no_None((None, 1))

    with pytest.raises(AssertionError):
        utils.check_all_or_no_None((1, None))


def test_read_value_or_warn():
    h = Header()
    h['a'] = 3

    assert utils.read_value_or_warn('a', h) == 3
    assert utils.read_value_or_warn(['b', 'a'], h) == 3
    with pytest.warns(UserWarning):
        utils.read_value_or_warn('c', h)


def test_WFC3WFC3_magnitude_zpt_reader():
    h = Header()
    h['TELESCOP'] = 'HST'
    h['INSTRUME'] = 'WFC3   '
    h['PHOTPLAM'] = 3
    h['PHOTFLAM'] = 2
    exp_result = -2.5 * np.log10(2) - 21.10 - 5 * np.log10(3) + 18.692

    assert np.isclose(utils.WFC3_magnitude_zpt_reader(h), exp_result)

    h['TELESCOP'] = 'WST'
    with pytest.raises(AssertionError):
        utils.WFC3_magnitude_zpt_reader(h)

    h['TELESCOP'] = 'HST'
    h['INSTRUME'] = 'WFC4   '
    with pytest.raises(AssertionError):
        utils.WFC3_magnitude_zpt_reader(h)

    h['INSTRUME'] = 'WFC3   '
    h.pop('PHOTPLAM')
    with pytest.raises(KeyError):
        utils.WFC3_magnitude_zpt_reader(h)


def test_read_zeropoint_magnitude():
    # first keyword
    h = Header()
    h['MAGZPT'] = 3
    assert utils.read_zeropoint_magnitude(h) == 3

    # second keyword
    h = Header()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', VerifyWarning)
        h['MAGZEROPOINT'] = 3
    assert utils.read_zeropoint_magnitude(h) == 3

    # use WFC3 fct
    h = Header()
    h['TELESCOP'] = 'HST'
    h['INSTRUME'] = 'WFC3   '
    h['PHOTPLAM'] = 3
    h['PHOTFLAM'] = 2
    exp_result = -2.5 * np.log10(2) - 21.10 - 5 * np.log10(3) + 18.692

    assert np.isclose(utils.WFC3_magnitude_zpt_reader(h), exp_result)

    # warn if not found
    h = Header()
    with pytest.warns(UserWarning):
        magzpt = utils.read_zeropoint_magnitude(h)
    assert magzpt is None
