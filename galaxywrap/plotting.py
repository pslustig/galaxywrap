import matplotlib.pyplot as plt
from .zscale import zscale


def plot_results(results, axes=None, *args, **kwargs):

    vmin, vmax = zscale(results['image'])
    kwargs['vmin'] = kwargs.pop('vmin', vmin)
    kwargs['vmax'] = kwargs.pop('vmax', vmax)
    kwargs['origin'] = 'lower'

    if axes is None:
        fig, axes = plt.subplots(1, 3, sharex=True, sharey=True)

    for ax, ext in zip(axes, ('image', 'model', 'residuals')):
        ax.imshow(results[ext], **kwargs)

    return axes
