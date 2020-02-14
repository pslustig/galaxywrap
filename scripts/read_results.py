import galaxywrap as gw
import argparse
from astropy.table import Table, Column, vstack
import numpy as np


def remove_columns_with_patterns(table, patterns):
    for k in table.keys():
        if any([p in k for p in patterns]):
            table.remove_column(k)


parser = argparse.ArgumentParser(description='print galfit results')
parser.add_argument('filenames', type=str, help='galfit imgblock filename(s)',
                    nargs='*', default=['imgblock.fits'])
parser.add_argument('-c', '--component', type=int, help='component to print')
parser.add_argument('-u', action='store_true', help='print uncertainties')
parser.add_argument('-f', action='store_true', help='print flags')
parser.add_argument('-s', action='store_true', help='print sky')
args = parser.parse_args()

patterns = [k for k, v in {
    '_unc': args.u, '_flag': args.f, 's': args.s}.items() if not v]

results = Table()
for i, filename in enumerate(args.filenames):
    result = gw.core.read_results(filename)['components']
    result.add_column(Column(np.arange(1, len(result)+1)), name='N', index=0)
    result.add_column(Column([i+1]*len(result)), name='fileno', index=0)
    result = Table(result[args.component-1]) if args.component else result
    results = vstack([results, result])

remove_columns_with_patterns(results, patterns)

results.pprint(max_lines=-1, max_width=-1, show_name=True,
               show_unit=None, show_dtype=False, align=None)
