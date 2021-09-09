import re
import warnings
from pyhdf.SD import SD
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pathlib import Path
from JSB_tools import mpl_hist, convolve_gauss, calc_background, convolve_gauss2d
from JSB_tools.TH1 import rolling_median
from JSB_tools.regression import PeakFit
from typing import Union, List, Dict, Tuple, Any
data_dir = Path(__file__).parent/'data'


def _load_array(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"file not found, {path}")

    f = SD(str(path))
    sds_obj = f.select('Streak_array')
    a = np.array(sds_obj, dtype=float)
    data0 = a[0]
    data1 = a[1]
    data = data0 - data1

    return data


def shot_metadata() -> Tuple[Dict[str, Any], List[str]]:
    with open(data_dir/'shot_metadata') as f:
        lines = f.readlines()
    lines = list(map(str.rstrip, lines))
    labels = {}
    out = {}
    for i in range(6):
        labels[i] = lines[i]
    _shot_num = 0
    for i in range(6, len(lines)):
        label_i = i % 6
        label = labels[label_i]
        text = lines[i]
        value = text
        if re.match('^[0-9]+$', value):
            value = int(value)
        elif re.match('^[0-9.]+$', value):
            value = float(value)
        if label == 'Shot':
            _shot_num = value + 35427 - 1
            out[_shot_num] = {label: _shot_num, 'True shot num':value}
        else:
            assert _shot_num != 0, 'Something went wrong'
            out[_shot_num][label] = text

    return out, list(labels.values())


# _shot_metadata, valid_labels = shot_metadata()
# for k, v in _shot_metadata.items():
#     print(k, v)


def get_shotnums_by_label(label):
    """
    Todo: Make this painfully complicated function simpler
    Args:
        label:

    Returns:

    """
    assert label in valid_labels, f'Invalid label, "{label}". Valid labels are {valid_labels}'
    m = {d[label]: d.copy() for _, d in _shot_metadata.items()}
    groups = {}
    for k, dat in m.items():
        full_dict = {k:v for k,v in dat.items()}
        dat.pop(label)
        dat.pop('Shot')
        dat.pop('True shot num')
        dat.pop('RID')
        key = tuple(dat.items())
        try:
            groups[key].append(full_dict)
        except KeyError:
            groups[key] = [full_dict]
    for key, list_of_infos in groups.items():
        list_of_infos = list(sorted(list_of_infos, key=lambda x: x[label]))
        groups[key] = list_of_infos
    out = {list_of_infos[0][label]: [l['Shot'] for l in list_of_infos] for list_of_infos in groups.values()}
    return out


class PJX2:
    def __init__(self, shot_number):
        self.shot_num = shot_number
        self.data = _load_array(data_dir / f'PJX-s{shot_number}_pjx2.hdf')
        self.data = np.sum(self.data, axis=0)  # spacial axis not interesting.. for now
        self._fit: Union[PeakFit, None] = None

    def align_peak(self, loc=1300):
        fit = self.fit_peak()
        center_index = np.searchsorted(np.arange(len(self.data)), float(fit.center.n))
        self.data = np.roll(self.data, shift=int(loc-center_index))

    @property
    def channels(self):
        return np.arange(len(self.data))

    def fit_peak(self,  convolve=10, bg_window=220, fit_window_width=100):
        if convolve:
            data = convolve_gauss(self.data, convolve)
        else:
            data = self.data

        if bg_window:
            data_bg = rolling_median(200, data)
            data -= data_bg

        peak_loc = np.argmax(data)
        _x = np.arange(len(data))
        peak_fit = PeakFit(peak_loc, _x, data, window_width=fit_window_width)
        self._fit = peak_fit
        return peak_fit

    @staticmethod
    def _get_label(_shot_num, labels: List[str]):
        strings = []
        try:
            meta_data = _shot_metadata[_shot_num]
        except KeyError:
            warnings.warn(f'No metadata for shot {_shot_num}. No label.')
            return ''
        for l in labels:
            assert l in valid_labels, f'Invalid label, "{l}". Valid labels are {valid_labels}'
            strings.append(f'{l}={meta_data[l]}')
        return '  '.join(strings)

    def plot_pjx2(self, ax=None, fig_size=6, convolve=10, bg_window=220, cut=(200, -200), labels=None,
                  plot_fit=False):
        """

        Args:
            fig_size:
            ax:
            convolve:
            bg_window:
            cut: The beginning and (maybe) end of the spectrum are not useful, make for ugly plots.
            labels:
            plot_fit:

        Returns:

        """
        if labels is None:
            labels = ['Shot', 'Pulse (ps)', 'Thick. (um)', 'Energy (J)']
        label = PJX2._get_label(self.shot_num, labels)
        if ax is None:
            plt.figure(figsize=(16/9*fig_size, fig_size))
            ax = plt.gca()

        if convolve:
            data = convolve_gauss(self.data, convolve)
        else:
            data = self.data

        if bg_window:
            data_bg = rolling_median(200, data)
            data -= data_bg
        data = data[slice(*cut)]

        peak_loc = np.argmax(data)
        _x = np.arange(len(data))
        peak = PeakFit(peak_loc, _x, data, window_width=100)
        ax, _ = mpl_hist(np.arange(len(data) + 1), data, label=f'{label}\namplitude: {peak.amp.n:.2e}  '
                                                               f'sigma: {peak.sigma: .2f}', ax=ax)
        if plot_fit:
            peak.plot_fit()

        ax.legend(prop={'size': 11})
        # ax.set_title(f"shot {self.shot_num}, PJX 2")
        ax.set_xlabel('Time channel')
        ax.set_ylabel('Intensity (integrated)')
        return ax


def plot_by_label(label):
    grouped_shots = get_shotnums_by_label(label)
    fig, axs = plt.subplots(len(grouped_shots))

    for ax, (value, shots) in zip(axs, grouped_shots.items()):
        for shot in shots:
            try:
                pjx = PJX2(shot)
            except FileNotFoundError:
                continue
            pjx.align_peak()
            pjx.plot_pjx2(ax=ax)


class CCD:
    def __init__(self, shot_num):
        self.shot_num = shot_num
        self.data = _load_array(data_dir/f'HRS_CCD-s{shot_num}_hrs_ccd_lower.hdf')
        self.roi = None

    def plot(self, ax=None):
        if ax is None:
            plt.figure()
            ax = plt.gca()

        data = convolve_gauss(np.sum(self.data[:, self.roi[0]: self.roi[1]], axis=1), 5)
        w = self.roi[1] - self.roi[0]
        s = 250
        # print(np.sum(data))
        data_bg1 = convolve_gauss(np.sum(self.data[:, self.roi[0]-w-s: self.roi[0]-s], axis=1), 5)
        data_bg2 = convolve_gauss(np.sum(self.data[:, self.roi[0] + s: self.roi[0] + w + s], axis=1), 5)
        print(list(map(lambda x: f'{np.sum(x):.2e}', [data_bg1, data_bg2, data])))
        mpl_hist(np.arange(len(data) + 1), data_bg1, ax=ax, label='bg1')
        mpl_hist(np.arange(len(data) + 1), data_bg2, ax=ax, label='bg2')
        mpl_hist(np.arange(len(data) + 1), data, ax=ax, label='sig')
        ax.legend()
        return ax

    def _find_line(self, n_sigma=1):
        a = convolve_gauss(np.sum(c.data, axis=0), 20)
        max_loc = np.argmax(a)
        fit = PeakFit(max_loc, np.arange(len(a)), a)
        self.roi = (max_loc-int(n_sigma*fit.sigma.n), max_loc+int(n_sigma*fit.sigma.n))
        return self.roi

    def imshow(self):
        plt.figure()
        ax = plt.gca()
        self._find_line()
        ax.imshow(self.data, norm=colors.LogNorm(), aspect='auto')
        ax.plot([self.roi[0]]*2, ax.get_ylim(), c='black')
        ax.plot([self.roi[1]]*2, ax.get_ylim(),  c='black')
        return ax

c = CCD(35437)

plt.plot(convolve_gauss(np.sum(c.data, axis=0), 8))

c.imshow()
c.plot()

# ==================================
# shotnums = [35430, 35432, 35434, 35436]
# n_overlays = 4
# ==================================
# print(get_shotnums_by_label('Thick. (um)'))
# plot_by_label('Thick. (um)')
#
# for index, shot_num in enumerate(shotnums):
#     if index%n_overlays == 0:
#         fig, ax = plt.subplots()
#     pjx = PJX2(shot_num)
#     pjx.align_peak()
#     pjx.plot_pjx2(ax=ax)

plt.show()