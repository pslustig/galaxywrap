from pathlib import Path
import logging
import subprocess
from astropy.io import fits
import os
import re
import numpy as np
from astropy.table import Table, vstack
from shutil import rmtree
import tempfile as tf
import sys
import warnings


def make_galfit_directory():
    newdir = tf.mkdtemp(prefix='galfit')
    return Path(newdir)


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
    directory = make_galfit_directory()

    make_galfit_files(feedme, image, psf, constraints, directory)

    cmd = [galfitcmd, 'galfit.feedme']

    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             cwd=directory, universal_newlines=True)

    print_galfit_output(popen, verbose)

    return_code = popen.wait()
    # if return_code:
    #     raise subprocess.CalledProcessError(return_code, cmd)


    # maybe just read everything if fit is done and just load mode if sources
    # are made
    results = read_results(directory)

    if deletefiles:
        rmtree(directory)

    if results is None:
        raise FitFailedError

    return results


def print_galfit_output(process, verbose):
    with process.stdout as pstd:
        for line in pstd:
            if verbose:
                print(line, end='')


class FitFailedError(Exception):
    def __init__(self, msg=None):
        if msg is None:
            msg = (r'''
                 __/~*##$%@@@******~\-__
               /f=r/~_-~ _-_ --_.^-~--\=b\
             4fF / */  .o  ._-__.__/~-. \*R\
            /fF./  . /- /' /|/|  \_  * *\ *\R\
           (iC.I+ '| - *-/00  |-  \  )  ) )|RB
           (I| (  [  / -|/^^\ |   )  /_/ | *)B
           (I(. \ `` \   \m_m_|~__/ )_ .-~ F/
            \b\\=_.\_b`-+-~x-_/ .. ,._/ , F/
             ~\_\= =  =-*###%#x==-#  *=- =/
                ~\**U/~  | i i | ~~~\===~
                        | I I \\
                       / // i\ \\
                  (   [ (( I@) )))  )
                       \_\_VYVU_/
                         || * |
                        /* /I\ *~~\
                      /~-/*  / \ \ ~~M~\
            ____----=~ // /WVW\* \|\ ***===--___

   Doh!  GALFIT crashed because at least one of the model parameters
   is bad.  The most common causes are: effective radius too small/big,
   component is too far outside of fitting region (also check fitting
   region), model mag too faint, axis ratio too small, Sersic index
   too small/big, Nuker powerlaw too small/big.  If frustrated or
   problem should persist, email for help or report problem to:
                     Chien.Y.Peng@gmail.com

        ''')
        super().__init__(msg)



def read_results(directory, filename='imgblock.fits'):
    if not check_if_fit_worked(directory):
        return None

    out = {}

    with fits.open(Path(directory)/filename) as hdul:
        header = fits.header.Header(hdul[2].header)
        out['image'] = image = hdul[1].data
        out['model'] = hdul[2].data
        out['residuals'] = hdul[3].data


    out['fitstats'] = read_fitstats_from_header(header)
    out['components'] = read_components_from_header(header)

    return out


def check_if_fit_worked(directory):
    imgblock_exists = (directory/'imgblock.fits').exists()
    return imgblock_exists

def read_fitstats_from_header(header):

    stats = {}
    stats['magzpt'] = header["MAGZPT"]

    fitreg = header["FITSECT"]
    fitreg = re.findall(r"[\w']+", fitreg)
    stats['box_x0'] = fitreg[0]
    stats['box_x1'] = fitreg[1]
    stats['box_y0'] = fitreg[2]
    stats['box_y1'] = fitreg[3]

    # Convolution box
    convbox = header["CONVBOX"]
    convbox = convbox.split(",")
    stats['convbox_x'] = convbox[0]
    stats['self.convbox_y'] = convbox[1]

    # Read in the chi-square value
    stats['chisq'] = header["CHISQ"]
    stats['ndof'] = header["NDOF"]
    stats['nfree'] = header["NFREE"]
    stats['reduced_chisq'] = header["CHI2NU"]

    return stats


def read_components_from_header(header):
    components_packs = make_component_packs(header)
    components = Table()
    for i, componentheader in enumerate(components_packs):
        components = vstack(
          [components, make_component_from_cleaned_header(componentheader, i)])

    return components


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
    translator = keywordtranslator()
    compname = header.pop('COMP_{}'.format(idx+1)).rstrip()

    t = Table()
    t['comp'] = [compname]
    for key, value in header.items():
        name = translator.to_python(key.replace('{}_'.format(idx+1), ''))
        add_parameter_to_table(t, name, *read_parameter(value))

    return t


def add_parameter_to_table(table, name, value, uncertainty, flag):
    ''' have to subtract 1 from coordinates due to indexing '''
    if name == 'x' or name == 'y':
        value -= 1
    table[name] = [value]
    table['{}_unc'.format(name)] = [uncertainty]
    table['{}_flag'.format(name)] = [flag]


def test_if_all_symbols_in_string(string, *symbols):
    isin = True
    for symbol in symbols:
        isin = isin and symbol in string

    return isin


def isconstrained(headerentry):
    return test_if_all_symbols_in_string(headerentry, '{', '}')


def isfixed(headerentry):
    return test_if_all_symbols_in_string(headerentry, '[', ']')


def isproblematic(headerentry):
    isproblematic = False
    if '*' in headerentry:
        isproblematic = True
    return isproblematic


def remove_and_split_string(string, *removechars):
    for char in removechars:
        string = string.replace(char, '')

    return string.split()


def read_parameter(headerentry):
    flag = ''
    uncertainty = -1

    if isconstrained(headerentry):
        flag = 'constrained'
        strval = remove_and_split_string(headerentry, '{', '}')

    elif isfixed(headerentry):
        flag = 'fixed'
        strval = remove_and_split_string(headerentry, '[', ']')

    elif isproblematic(headerentry):
        flag = 'problematic'
        strval = remove_and_split_string(headerentry, '*')
        uncertainty = np.nan
        warnings.warn('One parameter is problematic')

    else:
        strval = remove_and_split_string(headerentry)
        uncertainty = float(strval[2])

    value = float(strval[0])

    return value, uncertainty, flag


class keywordtranslator(object):
    def __init__(self):
        self.python = ['r', 'x', 'y']
        self.galfit = ['RE', 'XC', 'YC']

    @classmethod
    def to_python(cls, key, inverse=False):
        trans = cls()
        a = trans.galfit
        b = trans.python
        if inverse:
            a, b = b, a

        if key in a:
            key = b[a.index(key)]

        return key.lower()

    @classmethod
    def to_galfit(cls, key):
        return cls.to_python(key, inverse=True).upper()


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


def find_galfit_executable():
    galfit = ''
    try:
        galfit = os.environ['galfit']
    except KeyError:
        logging.warn(('Galfit not found in environment. Set path manually by '
                      'modifying the galfitcmd variable in core module.'))

    return galfit


galfitcmd = find_galfit_executable()
