import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import numpy as np

from . import utils


class WellPlotError(Exception):
    """
    Generic error class.
    """
    pass


def plot_depth_track_well(well, ax, md, kind='MD', tick_spacing=100):
    """
    Private function. Depth track plotting.

    Args:
        well (welly.well.Well): Well object.
        ax (ax): A matplotlib axis.
        md (ndarray): The measured depths of the track.
        kind (str): The kind of track to plot.

    Returns:
        ax.
    """
    if kind == 'MD':
        ax.set_yscale('bounded', vmin=md.min(), vmax=md.max())
    elif kind == 'TVD':
        tvd = well.location.md2tvd(md)
        ax.set_yscale('piecewise', x=tvd, y=md)
    else:
        raise Exception("Kind must be MD or TVD")

    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    for sp in ax.spines.values():
        sp.set_color('gray')

    if ax.is_first_col():
        pad = -10
        ax.spines['left'].set_color('none')
        ax.yaxis.set_ticks_position('right')
        for label in ax.get_yticklabels():
            label.set_horizontalalignment('right')
    elif ax.is_last_col():
        pad = -10
        ax.spines['right'].set_color('none')
        ax.yaxis.set_ticks_position('left')
        for label in ax.get_yticklabels():
            label.set_horizontalalignment('left')
    else:
        pad = -30
        for label in ax.get_yticklabels():
            label.set_horizontalalignment('center')

    ax.tick_params(axis='y', colors='gray', labelsize=12, pad=pad)
    ax.set_xticks([])

    ax.set(xticks=[])
    ax.depth_track = True

    return ax


def plot_well(well,
              legend=None,
              tracks=None,
              track_titles=None,
              alias=None,
              basis=None,
              return_fig=True,
              extents='td',
              **kwargs):
    """
    Plot multiple tracks.
    Args:
        well (welly.well.Well): Well object.
        legend (striplog.legend): A legend instance.
        tracks (list): A list of strings and/or lists of strings. The
            tracks you want to plot from ``data``. Optional, but you will
            usually want to give it.
        track_titles (list): Optional. A list of strings and/or lists of
            strings. The names to give the tracks, if you don't want welly
            to guess.
        alias (dict): a dictionary mapping mnemonics to lists of mnemonics.
        basis (ndarray): Optional. The basis of the plot, if you don't
            want welly to guess (probably the best idea).
        return_fig (bool): Whether to return the matplotlig figure. Default
            False.
        extents (str): What to use for the y limits:
            'td' — plot 0 to TD.
            'curves' — use a basis that accommodates all the curves.
            'all' — use a basis that accommodates everything.
            (tuple) — give the upper and lower explictly.

    Returns:
        None. The plot is a side-effect.
    """
    # These will be treated differently.
    depth_tracks = ['MD', 'TVD']

    # Set tracks to 'all' if it's None.
    tracks = tracks or list(well.data.keys())
    track_titles = track_titles or tracks

    # Check that there is at least one curve.
    if well.count_curves(tracks, alias=alias) == 0:
        if alias:
            a = " with alias dict applied "
        else:
            a = " "
        m = "Track list{}returned no curves.".format(a)
        raise WellPlotError(m)

    # Figure out limits
    if basis is None:
        basis = well.survey_basis(keys=tracks, alias=alias)

    if extents == 'curves':
        upper, lower = basis[0], basis[-1]
    elif extents == 'td':
        try:
            upper, lower = 0, well.location.td
        except:
            m = "Could not read well.location.td, try extents='curves'"
            raise WellPlotError(m)
        if not lower:
            lower = basis[-1]
    elif extents == 'all':
        raise NotImplementedError("You cannot do that yet.")
    else:
        try:
            upper, lower = extents
        except:
            upper, lower = basis[0], basis[-1]

    # Figure out widths because we can't us gs.update() for that.
    widths = [0.4 if t in depth_tracks else 1.0 for t in tracks]

    # Set up the figure.
    ntracks = len(tracks)
    fig = plt.figure(figsize=(2 * ntracks, 12), facecolor='w')
    fig.suptitle(well.name, size=16, zorder=100,
                 bbox=dict(facecolor='w', alpha=1.0, ec='none'))
    gs = mpl.gridspec.GridSpec(1, ntracks, width_ratios=widths)

    # Tick spacing
    order_of_mag = np.round(np.log10(lower - upper))
    ts = 10 ** order_of_mag / 100

    # Plot first axis.
    # kwargs = {}
    ax0 = fig.add_subplot(gs[0, 0])
    ax0.depth_track = False
    track = tracks[0]
    if '.' in track:
        track, kwargs['field'] = track.split('.')
    if track in depth_tracks:
        ax0 = well._plot_depth_track(ax=ax0, md=basis, kind=track, tick_spacing=ts)
    else:
        try:  # ...treating as a plottable object.
            ax0 = well.get_curve(track, alias=alias).plot(ax=ax0, legend=legend, **kwargs)
        except AttributeError:  # ...it's not there.
            pass
        except TypeError:  # ...it's a list.
            for t in track:
                try:
                    ax0 = well.get_curve(t, alias=alias).plot(ax=ax0, legend=legend, **kwargs)
                except AttributeError:  # ...it's not there.
                    pass
    tx = ax0.get_xticks()
    ax0.set_xticks(tx[1:-1])
    ax0.set_title(track_titles[0])

    # Plot remaining axes.
    for i, track in enumerate(tracks[1:]):
        # kwargs = {}
        ax = fig.add_subplot(gs[0, i + 1])
        ax.depth_track = False
        if track in depth_tracks:
            ax = well._plot_depth_track(ax=ax, md=basis, kind=track, tick_spacing=ts)
            continue
        if '.' in track:
            track, kwargs['field'] = track.split('.')
        plt.setp(ax.get_yticklabels(), visible=False)
        try:  # ...treating as a plottable object.
            curve = well.get_curve(track, alias=alias)
            curve._alias = track  # So that can retreive alias from legend too.
            ax = curve.plot(ax=ax, legend=legend, **kwargs)
        except AttributeError:  # ...it's not there.
            continue
        except TypeError:  # ...it's a list.
            for j, t in enumerate(track):
                if '.' in t:
                    track, kwargs['field'] = track.split('.')
                try:
                    curve = well.get_curve(t, alias=alias)
                    curve._alias = t
                    ax = curve.plot(ax=ax, legend=legend, **kwargs)
                except AttributeError:
                    continue
                except KeyError:
                    continue

        tx = ax.get_xticks()
        ax.set_xticks(tx[1:-1])
        ax.set_title(track_titles[i + 1])

    # Set sharing.
    axes = fig.get_axes()
    utils.sharey(axes)
    axes[0].set_ylim([lower, upper])

    # Adjust the grid.
    gs.update(wspace=0)

    # Adjust spines and ticks for non-depth tracks.
    for ax in axes:
        if not ax.depth_track:
            ax.set(yticks=[])
            ax.autoscale(False)
            ax.yaxis.set_ticks_position('none')
            ax.spines['top'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
            for sp in ax.spines.values():
                sp.set_color('gray')

    if return_fig:
        return fig
    else:
        return None
