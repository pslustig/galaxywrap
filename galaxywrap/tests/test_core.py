from galaxywrap import core
from astropy.io import fits
import pytest
import os
from distutils import dir_util
from pathlib import Path
import numpy as np
from astropy.table import Table


def get_directory():
    return Path(__file__).parent


def load_header():
    return fits.getheader(get_directory()/'test_core/imgblock.fits', 2)


def test_keywordtranslator():
    assert core.keywordtranslator.to_python('RE') == 'r'
    assert core.keywordtranslator.to_python('AR') == 'ar'
    assert core.keywordtranslator.to_python('MAG') == 'mag'

    assert core.keywordtranslator.to_galfit('r') == 'RE'
    assert core.keywordtranslator.to_galfit('ar') == 'AR'
    assert core.keywordtranslator.to_galfit('mag') == 'MAG'


def test_keywordtranslator_instance():
    trans = core.keywordtranslator()
    assert trans.to_python('RE') == 'r'
    assert trans.to_python('AR') == 'ar'
    assert trans.to_python('MAG') == 'mag'

    assert trans.to_galfit('r') == 'RE'
    assert trans.to_galfit('ar') == 'AR'
    assert trans.to_galfit('mag') == 'MAG'


def test_make_component_packs():
    h = load_header()
    p = core.make_component_packs(h)
    assert len(p) == 2

    assert len(p[0]) == 8
    assert len(p[1]) == 6

    assert p[0]['COMP_1'].rstrip() == 'sersic'
    assert '1_PA' in p[0].keys()

    assert p[1]['COMP_2'].rstrip() == 'sky'


def test_isconstrained():
    assert core.isconstrained('{46.5265 +/- 0.0220}')
    assert not core.isconstrained('46.5265 +/- 0.0220')


def test_isfixed():
    assert core.isfixed('[46.5265]')
    assert not core.isfixed('46.5265 +/- 0.0220')


def test_isproblematic():
    assert core.isproblematic('*46.5265* +/- *0.0220*')
    assert not core.isproblematic('46.5265 +/- 0.0220')


def test_remove_and_split_string():
    f = core.remove_and_split_string
    assert f('[]/a', '[', ']', '/') == ['a']
    assert f('[]/a') == ['[]/a']
    assert f('a b') == ['a', 'b']


def test_read_parameter():
    # why should error for constrained parameter be None?
    f = core.read_parameter
    assert f('{46.5265 +/- 0.0220}') == (46.5265, -1, 'constrained')
    assert f('[46.5265 +/- 0.0220]') == (46.5265, -1, 'fixed')
    r = f('*46.5265* +/- *0.0220*}')
    assert r[0::2] == (46.5265, 'problematic')
    assert np.isnan(r[1])


def test_add_parameter_to_table():
    t = Table()
    core.add_parameter_to_table(t, 'a', 1, 2, 'test')
    r = t[0]
    assert r['a'] == 1
    assert r['a_unc'] == 2
    assert r['a_flag'] == 'test'


def test_test_if_all_symbols_in_string():
    f = core.test_if_all_symbols_in_string
    assert f('{}', '{', '}')
    assert f('{}', '{')
    assert not f('{', '{', '}')
    assert not f('}', '{', '}')


def test_make_component_from_cleaned_header():
    h = load_header()
    p = core.make_component_packs(h)[0]
    core.make_component_from_cleaned_header(p, 0)


def test_read_components_from_header():
    h = load_header()
    print(core.read_components_from_header(h))

# test_make_component_from_cleaned_header()
test_read_components_from_header()
