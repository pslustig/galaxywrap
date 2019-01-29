from pathlib import Path
import logging
import subprocess
from astropy.io import fits
import os
import re


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
    deletefiles = kwargs.pop('deletefiles', True)
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


def read_results(directory, filename='imgblock.fits'):
    with fits.open(Path(directory)/filename) as hdul:
        header = fits.header.Header(hdul[2].header)
        model = hdul[2].data
        residuals = hdul[3].data

    ncomponents = get_number_of_component(header)

    return 0


def get_number_of_component(header):
    ncomps = 1
    while True:
        if "COMP_" + str(ncomps + 1) in header:
            ncomps = ncomps + 1
        else:
            break
    return ncomps


def make_component_packs(header):
    ''' takes whole header and appends all lines that belong so one single
        componend as a new item to a list '''
    components = []
    compidx = 1
    incomponent = False

    for i, (key, value) in enumerate(header.items()):
        _incomponent = is_part_of_component(key, compidx)

        if is_newcomponent(incomponent, _incomponent):
            startidx = i
        elif is_end_of_component(incomponent, _incomponent):
            components.append(header[startidx:i])
            compidx += 1

        incomponent = _incomponent

    return components


def make_component_from_cleaned_header(header, idx):
    parameters = {}
    uncertainties = {}
    modelname = header.pop('COMP_{}'.format(idx+1)).rstrip()
    for key, value in header.items():
        # read values and errors and put them in parameters and error dicts
        # 'flags?'
        pass


class keywordtranslator(object):
    def __init__(self):
        self.python = ['r', 'ratio']
        self.galfit = ['RE', 'RA']

    @classmethod
    def to_python(cls, key):
        trans = cls()
        if key in trans.galfit:
            key = trans.python[trans.galfit.index(key)]

        return key.lower()

    @classmethod
    def to_galfit(cls, key):
        trans = cls()
        if key in trans.python:
            key = trans.galfit[trans.python.index(key)]

        return key.upper()



def is_part_of_component(key, idx):
    ispart = False
    if key.startswith('{}_'.format(idx)) or (key == 'COMP_{}'.format(idx)):
        ispart = True
    return ispart


def is_newcomponent(incomponent_before, incomponent_now):
    isnew = False
    if (not incomponent_before) and incomponent_now:
        isnew = True

    return isnew


def is_end_of_component(incomponent_before, incomponent_now):
    isend = False
    if incomponent_before and (not incomponent_now):
        isend = True

    return isend

import re

class GalfitComponent(object):
    """Stores results from one component of the fit."""

    def __init__(self, galfitheader, component_number, verbose=True):
        """
        Read GALFIT results from output file.
        takes in the fits header from HDU 3 (the galfit model) from a
        galfit output file and the component number to extract
        """
        assert component_number > 0
        assert "COMP_" + str(component_number) in galfitheader

        self.component_type = galfitheader["COMP_" + str(component_number)]
        self.component_number = component_number
        headerkeys = [i for i in galfitheader.keys()]
        comp_params = []

        for i in headerkeys:
            if str(component_number) + '_' in i:
                comp_params.append(i)

        setattr(self, 'good', True)
        for param in comp_params:
            paramsplit = param.split('_')
            val = galfitheader[param]

            if "{" in val and "}" in val:
                print(" ## One parameter is constrained !")
                val = val.replace('{', '')
                val = val.replace('}', '')
                val = val.split()
                print(" ## Param - Value : ", param, val)
                setattr(self, paramsplit[1].lower(), float(val[0]))
                setattr(self, paramsplit[1].lower() + '_err', np.nan)
            elif "[" in val and "]" in val:
                print(" ## One parameter is fixed !")
                val = val.replace('[', '')
                val = val.replace(']', '')
                val = val.split()
                print(" ## Param - Value : ", param, val)
                setattr(self, paramsplit[1].lower(), float(val[0]))
                setattr(self, paramsplit[1].lower() + '_err', np.nan)
            elif "*" in val:
                print(" ## One parameter is problematic !")
                val = val.replace('*', '')
                val = val.split()
                print(" ## Param - Value : ", param, val)
                setattr(self, paramsplit[1].lower(), float(val[0]))
                setattr(self, paramsplit[1].lower() + '_err', -1.0)
                setattr(self, 'good', True)
            else:
                val = val.split()
                setattr(self, paramsplit[1].lower(), float(val[0]))
                setattr(self, paramsplit[1].lower() + '_err', float(val[2]))


class GalfitResults(object):

    """
    This class stores galfit results information.
    Currently only does one component
    """

    def __init__(self, galfitheader):
        """
        Init method for GalfitResults.
        Take in a string that is the name of the galfit output fits file
        """
        # Now some checks to make sure the file is what we are expecting
        # galfit_in_comments = False
        # for i in galfitheader['COMMENT']:
        #     galfit_in_comments = galfit_in_comments or "GALFIT" in i
        # assert True == galfit_in_comments
        # assert "COMP_1" in galfitheader
        # Now we've convinced ourselves that this is probably a galfit file

        # Read in the input parameters
        # self.input_initfile = galfitheader['INITFILE']
        # self.input_datain = galfitheader["DATAIN"]
        # self.input_sigma = galfitheader["SIGMA"]
        # self.input_psf = galfitheader["PSF"]
        # self.input_constrnt = galfitheader["CONSTRNT"]
        # self.input_mask = galfitheader["MASK"]
        # self.input_magzpt = galfitheader["MAGZPT"]

        # Fitting region
        # fitsect = galfitheader["FITSECT"]
        # fitsect = re.findall(r"[\w']+", fitsect)
        # self.box_x0 = fitsect[0]
        # self.box_x1 = fitsect[1]
        # self.box_y0 = fitsect[2]
        # self.box_y1 = fitsect[3]

        # Convolution box
        # convbox = galfitheader["CONVBOX"]
        # convbox = convbox.split(",")
        # self.convbox_x = convbox[0]
        # self.convbox_y = convbox[1]

        # Read in the chi-square value
        # self.chisq = galfitheader["CHISQ"]
        # self.ndof = galfitheader["NDOF"]
        # self.nfree = galfitheader["NFREE"]
        # self.reduced_chisq = galfitheader["CHI2NU"]
        # self.logfile = galfitheader["LOGFILE"]

        # Find the number of components
        num_components = 1
        while True:
            if "COMP_" + str(num_components + 1) in galfitheader:
                num_components = num_components + 1
            else:
                break
        self.num_components = num_components

        for i in range(1, self.num_components + 1):
            setattr(self, "component_" + str(i),
                    GalfitComponent(galfitheader, i),
                    )



galfitcmd = os.environ['galfit']
