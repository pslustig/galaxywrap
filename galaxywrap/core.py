from pathlib import Path
import logging
import subprocess
from astropy.io import fits
import os


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


def make_galfit_files(feedme, image, psf, constraints, directory):

    with open(directory / 'galfit.feedme', 'w') as outconf:
        outconf.write(feedme)

    fits.PrimaryHDU(image, header=image.properties.to_header()).writeto(
                                                        directory/'inimg.fits')
    if image.uncertainty is not None:
        fits.PrimaryHDU(image.uncertainty.array).writeto(
                                                directory/'sigma.fits')

    if psf is not None:
        fits.PrimaryHDU(psf.data).writeto(directory/'psf.fits')

    if image.mask is not None:
        fits.PrimaryHDU(image.mask.astype(int)).writeto(directory/'mask.fits')

    if constraints is not None:
        with open(directory/'constraints.txt', 'w') as outconst:
            outconst.write(constraints)

    return directory


def fit(feedme, image, psf, constraints, **kwargs):
    # add verbose
    verbose = kwargs.pop('verbose', False)
    directory = kwargs.pop('directory', '/tmp')
    directory = make_galfit_directory(directory)
    make_galfit_files(feedme, image, psf, constraints, directory)

    cmd = [galfitcmd, 'galfit.feedme']
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             universal_newlines=True, cwd=directory)
    '''
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    '''
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

    return directory


def read_results(directory):
    # galfit parser umbau, soll model einlesen
    return 0


galfitcmd = os.environ['galfit']
