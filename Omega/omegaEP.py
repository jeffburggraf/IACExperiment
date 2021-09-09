from pyhdf.SD import SD
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pathlib import Path
from JSB_tools import mpl_hist, convolve_gauss

"""
PJX 2 is Space and time
PJX 3 is energy and time
Material is buried to "tamp" it. The layer keeps target from expanding and cooling, thus makes temperature more uniform
SPC (Single-photon counting spectrometer) is just a normal camera. 
Look under the "DCHOPG" directory to get the dark room things.
CCD-time integrated, data is in an obvious stripe. 
"""
data_dir = Path(__file__).parent/'data'


def get_pjx_data(shot_num, pjx_num):
    if pjx_num == 2:
        fname = f'PJX-s{shot_num}_pjx2.hdf'
    elif pjx_num == 3:
        fname = f'PJX3-s{shot_num}_pjx3.hdf'
    else:
        assert False, 'Invalid PJX number'
    fpath = data_dir/fname
    if not fpath.exists():
        raise FileNotFoundError(f"File not found: '{fpath}' ")

    f = SD(str(fpath))
    # d = f.datasets()['Streak_array']
    sds_obj = f.select('Streak_array')  # select sds
    a = np.array(sds_obj, dtype=float)
    data0 = a[0]
    data1 = a[1]
    data = data0 - data1

    return data


def plot_pjx_data(shot_num, pjx_num, size=7):
    data = get_pjx_data(shot_num, pjx_num)
    plt.figure(figsize=(16/9*size, size))
    plt.imshow(data, norm=colors.LogNorm(), aspect='auto')
    plt.title(f"shot {shot_num}, PJX {pjx_num}")
    plt.xlabel('Time channel')
    plt.ylabel('Energy channel')
    plt.colorbar(label='intensity')

    return data, plt


streak_time_bins = {35429: (1476, 1502), 35430: (1460, 1490), 35431: (1486, 1514)}


def get_slice(shot_num, pjx_num, min_bin, max_bin, bg_offset=20, convolve=0):
    assert max_bin > min_bin, (min_bin, max_bin)
    w_width = max_bin - min_bin
    _d = get_pjx_data(shot_num, pjx_num=pjx_num)
    slice_sig = np.sum(_d[:, min_bin: max_bin], axis=1)/w_width
    slice_bg_left = np.sum(_d[:, min_bin-bg_offset-w_width: min_bin-bg_offset], axis=1)
    slice_bg_right = np.sum(_d[:, max_bin + bg_offset: max_bin + bg_offset + w_width], axis=1)
    slice_bg = (slice_bg_right + slice_bg_left) / (2*w_width)

    if convolve:
        slice_bg = convolve_gauss(slice_bg, convolve)
        slice_sig = convolve_gauss(slice_sig, convolve)
    # slice_sig -= slice_bg
    return slice_sig, slice_bg


def plot_slice(shot_num, pjx_num, min_bin=None, max_bin=None, bg_offset=20, bg=True, ax=None, convolve=0):
    if ax is None:
        plt.figure()
        ax = plt.gca()

    if min_bin is max_bin is None:
        min_bin, max_bin = streak_time_bins[shot_num]

    erg_slice_sig, erg_slice_bg = get_slice(shot_num, pjx_num, min_bin, max_bin, bg_offset, convolve=convolve)
    ax.set_title(f'From bin {min_bin} to {max_bin}')
    bins = np.arange(len(erg_slice_sig) + 1)
    mpl_hist(bins, erg_slice_sig, np.sqrt(erg_slice_sig), ax=ax, label=f"shot {shot_num}, PJX {pjx_num} (Signal)")
    if bg:
        mpl_hist(bins, erg_slice_bg, np.sqrt(erg_slice_bg), ax=ax, label=f"shot {shot_num}, PJX {pjx_num} (BG)")
    ax.legend()
    ax.set_xlabel("Energy bin")
    ax.set_ylabel("Intensity")
    return ax


#  ==================================
shot_nums = [35430, 35431]
pjx = 2
plot_bg = True
convolve = 5
#  ==================================
for shot_num in shot_nums:
    data, _ = plot_pjx_data(shot_num, 2, size=4)
    try:
        plot_slice(shot_num, 2, bg=plot_bg, convolve=convolve)
    except KeyError:
        print("No slice range, no slice")
plt.show()

