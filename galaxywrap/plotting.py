import matplotlib.pyplot as plt
from .zscale import zscale


def plot_results(results, *args, **kwargs):

    vmin, vmax = zscale(results['image'])
    kwargs['vmin'] = kwargs.pop('vmin', vmin)
    kwargs['vmax'] = kwargs.pop('vmax', vmax)
    kwargs['origin'] = 'lower'

    fig, axes = plt.subplots(1, 3, sharex=True, sharey=True)
    for i, ext in enumerate(('image', 'model', 'residuals')):
        axes[i].imshow(results[ext], **kwargs)

    return axes
