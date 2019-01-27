from pathlib import Path
import logging
import subprocess
from astropy.io import fits


def make_galfit_directory(where, exist_ok=False):
    where = Path(where)
    i = 0
    while True:
        subdir = where / 'galfit_{:03d}'.format(i)

        try:
            subdir.mkdir(parents=True, exist_ok=exist_ok)
            logging.info('Created directory %s', subdir)
            break
        except FileExistsError:
            i += 1

    return subdir


def make_directory_and_files(feedme, image, psf, mask, constraints, where):

    where = make_galfit_directory(where)

    with open(where / 'galfit.feedme', 'w') as outconf:
        outconf.write(feedme)

    fits.PrimaryHDU(image, header=image.to_header()).writeto(
                                                        where/'inimg.fits')
    if image.uncertainty is not None:
        fits.PrimaryHDU(image.uncertainty.array).writeto(where/'sigma.fits')

    if psf is not None:
        fits.PrimaryHDU(psf.data).writeto(where/'psf.fits')

    if image.mask is not None:
        fits.PrimaryHDU(image.mask.astype(int)).writeto(where/'mask.fits')

    if constraints is not None:
        with open(where/'constraints.txt', 'w') as outconst:
            outconst.write(constraints)

    return where


def fit(directory, verbose=False):
    cmd = ['galfit', 'galfit.feedme']
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def read_results():
    # galfit parser umbau, soll model einlesen
    pass
