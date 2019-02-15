import warnings
import numpy as np


def isiterable(obj):
    ''' apply iter on obj to see if it works and obj is iterable or not'''

    iterable = True
    try:
        iter(obj)
    except TypeError:
        iterable = False

    return iterable


def translate_to_constraints_names(name):
    if name == 'r':
        name = 're'
    return name


def check_all_or_no_None(tpl):
    if None in tpl:
        for entry in tpl:
            assert entry is None


def read_value_or_warn(keys, header):
    value = None

    if not isiterable(keys):
        keys = [keys]

    for key in keys:
        value = header.get(key, None)
        if value is not None:
            break

    if value is None:
        warnings.warn(
           'value with key(s) {} not found, must be set manually'.format(keys))

    return value


def WFC3_magnitude_zpt_reader(header):
    assert header['TELESCOP'] == 'HST'
    assert 'WFC3' in header['INSTRUME']
    photplam = header['PHOTPLAM']
    photflam = header['PHOTFLAM']
    return -2.5 * np.log10(photflam) - 21.10 - 5 * np.log10(photplam) + 18.692


def read_zeropoint_magnitude(header):
    keys = ['MAGZPT', 'MAGZEROPOINT']
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        magzpt = read_value_or_warn(keys, header)

    if magzpt is None:
        try:
            magzpt = WFC3_magnitude_zpt_reader(header)
        except (KeyError, AssertionError):
            warnings.warn('Could not find magnitude zeropoint, '
                          'must be set manually')

    return magzpt
