from analysis import Shot, SPEFile
import numpy as np
from matplotlib import pyplot as plt
from JSB_tools.MCNP_helper.outp_reader import StoppingPowerData
from JSB_tools import Nuclide, mpl_hist
from JSB_tools.list_reader import MaestroListFile


time_bins = np.arange(0, 50, 2)

attribs1 = {'flow': '100001',  'mylar': False,
            'foil_pos': 'upstream', 'flow_stop': None, 'num_filters': 2, 'beam_duration': 3}


def sum_shots(shot_nums=None, eval_func='1==1', remove_baseline=False, exclude_shots=None, **attribs):
    if exclude_shots is None:
        exclude_shots = []

    # used_shot_nums = []
    shot_nums = [] if shot_nums is None else shot_nums
    shots = [] if shot_nums is None else [Shot(i) for i in shot_nums]
    if len(attribs):
        for shot in Shot.find_shots(**attribs, eval_func=eval_func):
            if shot.shotnum in exclude_shots:
                continue
            if shot not in shots:
                shots.append(shot)

    l = None

    for shot in shots:
        print('Using shot: ', shot)
        if l is None:
            l = shot.list
        else:
            l += shot.list
    print()
    leg_label = f'Shots {",".join(map(str, shot_nums))}'

    l.path = leg_label
    return l, 1.0/len(shots)


# l1, s1 = sum_shots(shot_nums=[126, 127, 128])
# l2, s2 = sum_shots(shot_nums=[131, 132])
# MaestroListFile.multi_plotly([l1, l2], scales=[s1, s2], time_bin_width=6)

def warm_cold_filter_transit_time():
    """The transit times drop when filter is cold, all else equal. Why? """
    l1, s1 = sum_shots(shot_nums=[96, 96])  # Warm
    l2, s2 = sum_shots(shot_nums=[100, 101])    # cold
    MaestroListFile.multi_plotly([l1, l2], scales=[s1, s2],
                                 time_bin_width=6, erg_max=1600)



def mylar_tests():
    """
    Plot spectra from mylar ranging from 0 to 20 um.
    Returns:

    """
    # two shots for each My thickness and flow rate. Each shot is for different flow congif: 010010 or 100001
    shots_25 = {0: [40, 122], 2.5: [103, 105], 5: [107, 108], 10: [111, 133, 114], 20: [116, 118]}  # 0.25 L/m
    shots_5 = {0: [38, 121], 2.5: [102, 104], 5: [106, 109], 10: [110, 112], 20: [115, 117]}  # 0.5 L/m

    for d, flow in [(shots_25, '0.25'), (shots_5, '0.5')]:
        scales = []
        ls = []
        labels = []
        for k, v in shots_25.items():
            l, s = sum_shots(v)
            ls.append(l)
            scales.append(s)
            labels.append(f"{k} um My ({flow} L/m)")
        MaestroListFile.multi_plotly(ls, scales=scales,
                                     leg_labels=labels,
                                     time_bin_width=6, erg_max=1600)


def fresh_filters(fresh_only=False):
    """
    Only plot lines at hte beginning of the day (fresh filter).
    Shot 93 is a fresh filter.
    Should I include shot 136?
    Returns:

    """
    fresh = [1, 93, 124]
    day_beginning = [120, 83, 45]

    shots = fresh

    if not fresh_only:
        shots += day_beginning

    l = None
    for shot in map(Shot, shots):
        if l is None:
            l = shot.list
        else:
            l += shot.list

    l.plotly(erg_max=1800, remove_baseline=True, convolve_overlay_sigma=1)


#  Todo: flow pattern

# mylar_tests()
# warm_cold_filter_transit_time()
fresh_filters()
