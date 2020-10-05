import warnings
import numpy as np
import astropy.units as u
from scipy import signal
from astropy.convolution import Tophat2DKernel
from astropy import wcs


def cut_to_fitarea(array, fitarea):
    (x0, x1), (y0, y1) = fitarea
    return array[y0:y1+1, x0:x1+1]


def random_unmasked_positions(mask, npos, replace=True):
    ''' return npos pixel idices of randomly choosen unmasked pixels in order
        x, y '''
    ny, nx = mask.shape
    yx = np.array(np.meshgrid(np.arange(ny), np.arange(nx))).T
    unmasked = yx[~mask]
    randoms = unmasked[np.random.choice(len(unmasked), npos, replace=replace)]
    return randoms[:, 1], randoms[:, 0]


def img_to_numpy(img, masked=False, addmask=None):
    img = np.ma.masked_array(img.data, mask=img.mask)

    if addmask is not None:
        img.mask = img.mask | addmask

    if not masked:
        return img.filled(np.inf)
    return img


def grow_mask(mask, npix):
    ''' extend mask by npix '''
    kernel = Tophat2DKernel(npix)
    kernel.normalize('peak')
    return signal.fftconvolve(mask, kernel, mode='same') > 0.5


def shrink_mask(mask, npix):
    return ~grow_mask(~mask, npix)


def make_sourcemask(sourcecoords, WCS, maskradius):
    mask = np.zeros(WCS.array_shape, dtype=bool)
    mask += cat_on_grid(sourcecoords, WCS).astype(bool)

    pixscale = get_pixelscale(WCS)

    return grow_mask(mask, maskradius.to_value(u.pix, equivalencies=pixscale))


def get_pixelscale(WCS):
    return u.pixel_scale(np.mean(
            wcs.utils.proj_plane_pixel_scales(WCS)) * u.deg / u.pix)


def cat_on_grid(sourcecoords, gridwcs):
    radecbin = sourcecoords.to_pixel(gridwcs)
    hist, _, _ = np.histogram2d(*radecbin, bins=[
                                np.arange(gridwcs.array_shape[1]+1),
                                np.arange(gridwcs.array_shape[0]+1)])
    return hist.T


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
