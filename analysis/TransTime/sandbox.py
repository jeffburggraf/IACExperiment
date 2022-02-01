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

l1, s1 = sum_shots(shot_nums=[96, 96])
l2, s2 = sum_shots(shot_nums=[100, 101])
MaestroListFile.multi_plotly([l1, l2], scales=[s1, s2], time_bin_width=6, erg_max=1600)



