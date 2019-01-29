from pathlib import Path
import galaxywrap as gw
from astropy.io import fits
import warnings
from importlib import reload
from astropy.nddata import StdDevUncertainty

reload(gw)

%load_ext autoreload
%autoreload 2


datadir = Path().home() / 'Documents/galaxywrap/galaxywrap/tests/files'

with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)
    img = gw.image.read(datadir/'galaxy.fits')

img.uncertainty = StdDevUncertainty(fits.getdata(datadir/'uncertainty.fits'))

psf = gw.psf(fits.getdata(datadir/'psf.fits'), (100, 100))
component = gw.models.sersic(48.518, 51.2800, 20.0890, 5.1160, 4.2490, 0.7570, -60.3690)
model = gw.models.model(component)

print(model.fit(img, psf, ''))
