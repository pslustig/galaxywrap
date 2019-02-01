from pathlib import Path
import galaxywrap as gw
from astropy.io import fits
import warnings


datadir = Path().home() / 'Documents/galaxywrap/galaxywrap/tests/files' # mac
datadir = Path().home() / 'projects/galaxywrap/galaxywrap/tests/files' # uni

with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)
    img = gw.image.read(datadir/'gal.fits')

psf = gw.psf(fits.getdata(datadir/'psf.fits'), (100, 100))
component = gw.models.sersic(48.518, 51.2800, 20.0890, 5.1160, 4.2490, 0.7570, -60.3690)
model = gw.models.model(component)
model.add_component(gw.models.sky(1.3920, 0, 0,
                                  fixed={'dbkg_dx': True, 'dbkg_dy': True}))

res = model.fit(img, psf, '')
