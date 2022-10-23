from analysis import Shot, SPEFile
import numpy as np
from matplotlib import pyplot as plt
from JSB_tools.MCNP_helper.outp_reader import StoppingPowerData
from JSB_tools import Nuclide, mpl_hist
from JSB_tools.maestro_reader import MaestroListFile


time_bins = np.arange(0, 50, 2)

attribs1 = {'flow': '100001',  'mylar': False,
            'foil_pos': 'upstream', 'flow_stop': None, 'num_filters': 2, 'beam_duration': 3}


def sum_shots(shot_nums=None, eval_func='1==1', remove_baseline=False, exclude_shots=None, **attribs):
    if remove_baseline:
        raise NotImplementedError

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
            l.__iadd__(shot.list, truncate_time=True)
    print()
    leg_label = f'Shots {",".join(map(str, shot_nums))}'

    # l.path = leg_label
    return l, 1.0/len(shots), leg_label


# l1, s1 = sum_shots(shot_nums=[126, 127, 128])
# l2, s2 = sum_shots(shot_nums=[131, 132])
# MaestroListFile.multi_plotly([l1, l2], scales=[s1, s2], time_bin_width=6)

def warm_cold_filter_transit_time(time_bin_width=6, he_flow=0.5):
    """The transit times drop when filter is cold, all else equal. Why?

    Consider making a time dependence plot showing this effect for paper.

     """
    # if he_flow == 0.5:
    #     warm_shots = [93, 94]
    #     cold_shots = [98, 99]
    # elif he_flow == 1.0:
    #     assert False
    # else:
    #     assert False
    warm_shots = [124, 125, 126, 127]
    cold_shots = [129, 130, 131, 132]
    l1, s1, leg_label1 = sum_shots(shot_nums=cold_shots)    # cold
    l2, s2, leg_label2 = sum_shots(shot_nums=warm_shots)  # Warm
    MaestroListFile.multi_plotly([l1, l2], scales=[s1, s2],
                                 time_bin_width=time_bin_width, erg_max=1600,
                                 leg_labels=[f"Cold-{leg_label1}", f"Warm-{leg_label2}"])


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
            l, s, _ = sum_shots(v)
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
        print(shot)
        if l is None:
            l = shot.list
        else:
            l += shot.list

    l.plotly(erg_max=1800, remove_baseline=True, convolve_overlay_sigma=1)


def all_fast_warm_shots(**plotly_kwargs):
    """
    All warm filter shots with short tub and 0.5/0.5 Ar/He flow rate.
    Returns:

    """
    list_ = None

    for shot in Shot.find_shots(tube_len=4.16, cold_filter=False, mylar=0, flow_stop=0,  beam_duration=3,
                                num_filters=2, eval_func='self.he_flow + self.ar_flow >= 1.0', ):
        print(shot)
        l = shot.list
        if list_ is None:
            list_ = l
        else:
            list_ += l

    list_.plotly(**plotly_kwargs, time_bin_width=7, time_step=3)


def cold_v_warm():
    cold_shots = range(129, 135)
    warm_shots = range(124, 129)
    warm_list = None
    cold_list = None

    for s in cold_shots:
        if cold_list is None:
            cold_list = Shot(s).list
        else:
            cold_list += Shot(s).list

    for s in warm_shots:
        if warm_list is None:
            warm_list = Shot(s).list
        else:
            warm_list += Shot(s).list

    MaestroListFile.multi_plotly([cold_list, warm_list], leg_labels=['Cold', 'Warm'])


def inline_vs_offset(flow_rate=0.5, tube_len=4.16, max_shots=5, foil_pos='center', baseline_subtract=False):
    offset_list = None
    inline_list = None

    n_inline = 0
    n_offset = 0

    for shot in Shot.find_shots(tube_len=tube_len, cold_filter=False, mylar=0, flow_stop=0, beam_duration=3,
                                num_filters=2, he_flow=flow_rate, ar_flow=flow_rate, flow='010010', foil_pos=foil_pos):
        print('Inline:', shot)
        if n_inline > max_shots:
            break

        n_inline += 1
        if inline_list is None:
            inline_list = shot.list
        else:
            inline_list.__iadd__(shot.list, truncate_time=True)
        # inline_list += shot.list

    for shot in Shot.find_shots(tube_len=tube_len, cold_filter=False, mylar=0, flow_stop=0, beam_duration=3,
                                num_filters=2, he_flow=flow_rate, ar_flow=flow_rate, flow='100001', foil_pos=foil_pos):
        print('Offset:', shot)
        n_offset += 1
        if offset_list is None:
            offset_list = shot.list
        else:
            offset_list.__iadd__(shot.list, truncate_time=True)

        if n_offset > max_shots:
            break

    MaestroListFile.multi_plotly([offset_list, inline_list],
                                 leg_labels=['Offset', 'Inline'],
                                 scales=[1.0/n_offset, 1.0/n_inline],
                                 remove_baseline=baseline_subtract)

# Todo: flow pattern

inline_vs_offset()
# cold_v_warm()

# mylar_tests()
# all_fast_warm_shots()
# warm_cold_filter_transit_time(15)

# fresh_filters()
