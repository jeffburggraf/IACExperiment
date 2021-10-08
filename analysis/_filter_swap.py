import matplotlib.pyplot as plt
import numpy as np
from analysis import Shot, get_merged_time_dependence
# import numpy as np
from JSB_tools import mpl_hist, convolve_gauss
# import matplotlib.pyplot as plt
from JSB_tools.regression import LinearFit
from uncertainties import ufloat
from uncertainties import unumpy as unp
from lmfit.models import LinearModel
from lmfit.models import ExponentialModel


def filter_consistency():
    shots1 = [79, 85, 86]
    shots2 = [83, 84]
    ax = None
    sigs = []
    for shotnum in shots1:
        shot = Shot(shotnum)
        y, _, bins = shot.list.get_time_dependence(218.5)
        sigs.append(sum(y))
        sigs[-1] = ufloat(sigs[-1], np.sqrt(sigs[-1]))
        if ax is None:
            ax = mpl_hist(bins, y, color='blue')
        else:
            mpl_hist(bins, y, color='blue', ax=ax)

    for shotnum in shots2:
        shot = Shot(shotnum)
        y, _, bins = shot.list.get_time_dependence(218.5)
        sigs.append(sum(y))
        sigs[-1] = ufloat(sigs[-1], np.sqrt(sigs[-1]))
        mpl_hist(bins, y, color='green', ax=ax)

    print(f"Standard Deviation between shots:{100*np.std(unp.nominal_values(sigs))/np.mean(unp.nominal_values(sigs)):.1f}%"
          f" Counting errors: {100*np.mean(unp.std_devs(sigs))/np.mean(unp.nominal_values(sigs)):.1f}%")


def hi(center_kev=218.6, window_kev=3, bins=None, max_time=None, debug_plot=False):
    shots_upstream = [79, 85, 86]
    shots_downstream = [83, 84]
    # Shot(132).list.plotly(100, 900, time_bin_width=40)

    upstream = None
    downstream = None
    erg_bins = None
    erg_downstream = None
    erg_upstream = None
    down_time_tot = 0
    up_time_tot = 0
    min_time_up = None # min time of upstream shots
    min_time_down = None  # min time of downstream shots

    for shot_up in shots_upstream:
        shot = Shot(shot_up)
        print(f'shot {shot_up} livetime: {shot.list.total_livetime}')

        spe = shot.list.build_spe(max_time=max_time)
        up_time_tot += shot.list.total_livetime
        if min_time_up is None or shot.list.total_livetime < min_time_up:
            min_time_up = shot.list.total_livetime
        if bins is None:
            _upstream, _, bins = shot.list.get_time_dependence(center_kev, signal_window_kev=window_kev,
                                                               debug_plot=debug_plot)
        else:
            _upstream, _, _ = shot.list.get_time_dependence(center_kev, signal_window_kev=window_kev, bins=bins,
                                                            debug_plot=debug_plot)

        if erg_upstream is None:
            erg_bins = spe.erg_bins
            erg_upstream = spe.get_counts(remove_baseline=True)
        else:
            erg_upstream += spe.get_counts(remove_baseline=True)

        if upstream is None:
            upstream = _upstream
        else:
            upstream += _upstream

    for shot_down in shots_downstream:
        shot = Shot(shot_down)
        print(f'shot {shot_down} livetime: {shot.list.total_livetime}')

        down_time_tot += shot.list.total_livetime
        spe = shot.list.build_spe(max_time=max_time)

        if min_time_down is None or shot.list.total_livetime < min_time_down:
            min_time_down = shot.list.total_livetime

        _downstream, _, _ = shot.list.get_time_dependence(center_kev, signal_window_kev=window_kev, bins=bins,
                                                          debug_plot=debug_plot)
        if downstream is None:
            downstream = _downstream
        else:
            downstream += _downstream

        if erg_downstream is None:
            erg_downstream = spe.get_counts(remove_baseline=True)
        else:
            erg_downstream += spe.get_counts(remove_baseline=True)

    ax = mpl_hist(bins, downstream*sum(upstream)/sum(downstream), label=f'downstream (scaled), shots:{shots_downstream}')
    mpl_hist(bins, upstream, label=f'upstream shots:{shots_upstream}', ax=ax)
    ax.set_title(f"Comparison of time dependence in first/seconds filters. {center_kev} KeV")

    ax = mpl_hist(bins, downstream, label='downstream (not normalized)')
    mpl_hist(bins, upstream, label='upstream', ax=ax)
    ax.set_title(f"Comparison of rates in first/seconds filters. {center_kev} KeV")

    ax = mpl_hist(erg_bins, erg_upstream/up_time_tot, label=f'upstream ({up_time_tot:.0f}s)')
    mpl_hist(erg_bins, erg_downstream/down_time_tot, label=f'down stream ({down_time_tot:.0f}s)', ax=ax)
    bg = Shot.background_spe()
    mpl_hist(bg.erg_bins, bg.get_counts(remove_baseline=True, make_rate=True), ax=ax, label='bg')

    print(f'Min_time_down: {min_time_down}; min_time_up: {min_time_up}')
    print(f"Ratio: {sum(downstream)/ sum(upstream)}")


bins = np.arange(0, 450, 7)
hi(center_kev=743, bins=bins, max_time=bins[-1], debug_plot='simple')
hi(center_kev=229, bins=bins, max_time=bins[-1], debug_plot='simple')
hi(center_kev=658, bins=bins, max_time=bins[-1], debug_plot='simple')



plt.show()

