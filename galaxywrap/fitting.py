from . import utils
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.table import Table, Row
from . import galaxywrap
from . import components


'''def prepare_fit(coord, WCS, size, sources, modelmaker=make_fitmodel_with_stars,
                *args, **kwargs):
    xc, yc = coord.to_pixel(WCS)
    calculate_image_coordinates(sources, WCS)
    pixscale = utils.get_pixelscale(WCS)
    size_pix = size.to_value(u.pix, equivalencies=pixscale)

    fitzone = get_fitzone((xc, yc), size_pix, WCS.array_shape)

    sources = select_targets_in_cutout(fitzone, sources, buffer=5)
    return modelmaker(sources), fitzone
'''


def add_sersic_to_img(img, parameters, addnoise=True):
    parnames = 'x', 'y', 'mag', 'r', 'n', 'ar', 'pa'
    component = components.sersic(**{p: parameters[p] for p in parnames})
    return img.add_component(component)


def make_component_from_parameters(params):
    ''' dict or astropy.table.Row containing modelparameters and
        key comp containing the component name '''
    if isinstance(params, Row):
        params = {k: params[k] for k in params.colnames}
    return getattr(components, params.pop('comp'))(**params)


def make_fitmodel_from_sourcetable(sourcetable, fitarea=None):
    if fitarea is not None:
        sourcetable = select_targets_in_fitarea(fitarea, sourcetable)

    mdl = galaxywrap.model()
    for source in sourcetable:
        mdl.add_component(make_component_from_parameters(source))
    return mdl


def sextractor_to_galfit_table(t, WCS=None):
    x, y = t['X_IMAGE'], t['Y_IMAGE']
    mag, n, r = t['MAG_AUTO'], [1.5]*len(t), t['FLUX_RADIUS']
    ar = t['B_IMAGE'] / t['A_IMAGE']
    pa = t['THETA_IMAGE'] + 90
    comp = ['sersic'] * len(t)

    if WCS is not None:
        x, y = SkyCoord(t['X_WORLD'], t['Y_WORLD'], unit='deg').to_pixel(WCS)

    gf = Table(names=('x', 'y', 'mag', 'n', 'r', 'ar', 'pa', 'comp'),
               data=(x, y, mag, n, r, ar, pa, comp))

    if 'star' in t.keys():
        gf['comp'][t['star']] = 'psf'

    return gf


def make_artificial_source_table(sourceparameter, notcleanmask, nsources):
    positions = np.array(utils.random_unmasked_positions(
                        notcleanmask, nsources), dtype=float)
    positions += np.random.uniform(low=-.5, high=.5, size=positions.shape)
    pa = np.random.uniform(0, 90, nsources)
    poskeys = ('x', 'y', 'pa')
    t = Table(names=poskeys, data=(positions[0], positions[1], pa))
    for key, value in sourceparameter.items():
        if key not in poskeys:
            t[key] = value
    return t


def get_fitarea(center, fitsize, imgshape=None, WCS=None):
    ''' fitsize must be integer, same size in x and y is assumed '''
    pixscale = None if WCS is None else utils.get_pixelscale(WCS)
    if isinstance(center, SkyCoord):
        center = center.to_pixel(WCS)
    fitsize = fitsize.to_value(u.pix, equivalencies=pixscale)

    # BUG: fitzone does not result in requested shape (1 pixel too large)
    # either error here or in galfit fitzone wrapper
    xcenter, ycenter = center
    halfsize = fitsize / 2
    xmin, xmax = int(xcenter - halfsize), int(xcenter + halfsize)
    ymin, ymax = int(ycenter - halfsize), int(ycenter + halfsize)

    # crop to shape
    if imgshape is None:
        imgshape = WCS.array_shape

    xmin = max(xmin, 0)
    xmax = min(xmax, imgshape[1]-1)
    ymin = max(ymin, 0)
    ymax = min(ymax, imgshape[0]-1)

    return ((xmin, xmax), (ymin, ymax))


def select_targets_in_fitarea(fitarea, sourcetable, buffer=0):
    (x0, x1), (y0, y1) = fitarea
    x, y = sourcetable['x'], sourcetable['y']
    take = ((x > (x0 + buffer)) & (x < (x1 - buffer))
            & (y > (y0 + buffer)) & (y < (y1 - buffer)))
    return sourcetable[take]
