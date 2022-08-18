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
from matplotlib import rc

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
    for i in range(7):
        labels[i] = lines[i]
    _shot_num = 0
    for i in range(7, len(lines)):
        label_i = i % 7
        label = labels[label_i]
        text = lines[i]
        value = text
        if label == 'Target':
            value = value.lstrip('-')
        if re.match('^[0-9]+$', value):
            value = int(value)
        elif re.match('^[0-9.]+$', value):
            value = float(value)
        if label == 'Shot':
            _shot_num = value
            out[_shot_num] = {label: _shot_num}
        else:
            assert _shot_num != 0, 'Something went wrong'
            out[_shot_num][label] = value

    return out, list(labels.values())


_shot_metadata, valid_labels = shot_metadata()
for k, v in _shot_metadata.items():
    print(k, v)


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
    def _get_label(_shot_num, labels: List[str], bold_labels=None):
        if bold_labels is None:
            bold_labels = []
        else:
            rc('text', usetex=True)
        strings = []
        try:
            meta_data = _shot_metadata[_shot_num]
        except KeyError:
            warnings.warn(f'No metadata for shot {_shot_num}. No label.')
            return ''
        for l in labels:
            assert l in valid_labels, f'Invalid label, "{l}". Valid labels are {valid_labels}'
            if meta_data[l] is bold_labels:
                strings.append(f'{l}=' + r'\textbf{{{}}}'.format(meta_data[l]))
            else:
                strings.append(f'{l}={meta_data[l]}')
        return '  '.join(strings)

    def plot_pjx2(self, ax=None, fig_size=6, convolve=10,
                  bg_window=220, cut=(200, -200),
                  label: Union[str, List[str], None] = None,
                  plot_fit=False, fit_label=False,
                  bold_labels: List[str] = None, title=None, **plt_args):
        """

        Args:
            fig_size:
            ax:
            convolve:
            bg_window:
            cut: The beginning and (maybe) end of the spectrum are not useful, make for ugly plots.
            label: If None, use all shot config.
                   If list of strings, use only shot configs corresponding to the strings.
                    Valid labels: ['Shot', 'RID', 'Target', 'SL/BL', 'Pulse (ps)', 'Energy (J)', 'Thick. (um)']
                   If string, simply use that string.
                   If "None", then no label.
            plot_fit:
            fit_label: If True, add peak amplitude and sigma to label.
            bold_labels: Labels that wil be bold text.
            title:
            plt_args: Arguments to be passed to ax.plot()


        Returns:

        """
        if not isinstance(label, str):
            if label is None:
                labels = ['Shot', 'Pulse (ps)', 'Thick. (um)', 'Energy (J)']
            else:
                assert all(map(lambda x: isinstance(x, str), label)), 'Invalid `label` argument'
                labels = label
            label = PJX2._get_label(self.shot_num, labels, bold_labels=bold_labels)
        else:
            assert isinstance(label, str), 'Invalid label type. Must be of type Union[str, List[str], None]'
            if label == 'None':
                label = None
        if plt_args is None:
            plt_args = {}

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
        if fit_label:
            label += f'\namplitude: {peak.amp.n:.2e} sigma: {peak.sigma: .2f}'
        ax, _ = mpl_hist(np.arange(len(data) + 1), data, label=label, ax=ax, **plt_args)
        if plot_fit:
            peak.plot_fit()

        ax.legend(prop={'size': 11})
        ax.set_xlabel('Time channel')
        ax.set_ylabel('Intensity')
        if title is None:
            ax.set_title(f'PJ-X 2')
        else:
            ax.set_title(title)
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

    def plot(self, ax=None, convolve=35, label=None, title=None, **plt_args):
        if ax is None:
            plt.figure()
            ax = plt.gca()
        ax.set_title(f'CCD (ROI)-{self.shot_num}')
        data = np.sum(self.data[:, self.roi[0]: self.roi[1]], axis=1)
        data = rolling_median(convolve, data)
        # data = convolve_gauss(, convolve)
        w = self.roi[1] - self.roi[0]
        s = 250
        # print(np.sum(data))
        # data_bg1 = np.sum(self.data[:, self.roi[0] - w - s: self.roi[0] - s], axis=1)
        # data_bg1 = convolve_gauss(data_bg1, convolve)
        # data_bg2 = np.sum(self.data[:, self.roi[0] + s: self.roi[0] + w + s], axis=1)
        # data_bg2 = convolve_gauss(data_bg2, convolve)
        # data_bg1 = rolling_median(convolve, data_bg1)
        # data_bg2 = rolling_median(convolve, data_bg2)
        # mpl_hist(np.arange(len(data) + 1), data_bg1, ax=ax, label='bg1')
        # mpl_hist(np.arange(len(data) + 1), data_bg2, ax=ax, label='bg2')
        mpl_hist(np.arange(len(data) + 1), data, ax=ax, label=label, title=title, **plt_args)
        ax.legend()
        return ax
    
    def _find_line(self, n_sigma=1):
        a = convolve_gauss(np.sum(self.data, axis=0), 20)
        max_loc = np.argmax(a)
        fit = PeakFit(max_loc, np.arange(len(a)), a, fix_center=True)
        # fit.plot_fit()
        self.roi = (max_loc-int(n_sigma*fit.sigma.n), max_loc+int(n_sigma*fit.sigma.n))
        return self.roi

    def imshow(self):
        plt.figure()
        plt.title(f'CCD-{self.shot_num}')
        ax = plt.gca()
        self._find_line()
        ax.imshow(self.data, norm=colors.LogNorm(), aspect='auto')
        ax.plot([self.roi[0]]*2, ax.get_ylim(), c='black')
        ax.plot([self.roi[1]]*2, ax.get_ylim(),  c='black')
        return ax


def get_same_config():
    shot_dicts = {k: v.copy() for k, v in _shot_metadata.items()}
    outs = {}
    configs = []

    for shot_primary, d in shot_dicts.items():
        d.pop('Shot')
        d.pop('RID')
        d.pop('Target')
        key = tuple(d.items())
        for index, config in enumerate(configs):
            if d == config:
                outs[key].append(shot_primary)
                break
        else:
            outs[key] = [shot_primary]
            configs.append(d)
    return outs


print("Same configs:")
for k, v in get_same_config().items():
    print(k, v)


# ==================================
shotnums = [35427, 35428, 35429, 35431]
n_overlays = 4
# ==================================
fig, ax = plt.subplots()

labels = ['Energy (J) 250', 'Energy (J) 250', 'Energy (J) 500', 'Energy (J) 900']

for index, shot_num in enumerate(shotnums):
    if index%n_overlays == 0:
        fig, ax = plt.subplots()
#     pjx = PJX2(shot_num)
#     pjx.align_peak()
#     pjx.plot_pjx2(ax=ax)
    c = CCD(shot_num)
    c.imshow()
    ax.set_xlabel('Energy channel')
    ax.set_ylabel('Intensity')
    if index == 0:
        c.plot(ax, label=labels[index], title='CCD', color='red')
    elif index == 1:
        c.plot(ax, label=labels[index], title='CCD', color='red', ls='--')
    else:
        c.plot(ax, label=labels[index], title='CCD')


def temp(shot1, shot2, color, title):
    if shot1 is not None:
        pjx = PJX2(shot1)
        pjx.align_peak()
        pjx.plot_pjx2(ax=ax, label=labels, c=color, title=title)
    if shot2 is not None:
        pjx = PJX2(shot2)
        pjx.align_peak()
        pjx.plot_pjx2(ax=ax, label='None', c=color, title=title)

# labels = ['Thick. (um)']

#  ('Pulse (ps)', 'BC'), ('Energy (J)', 500), ('Thick. (um)', varies))
# temp(35430, None, 'blue', title='PJ-x 2; BL; Energy (J): 500 J')
# temp(35432, 35439, 'red', title='PJ-x 2; BL; Energy (J): 500 J')
# temp(35434, 35436, 'green', title='PJ-x 2; BL; Energy (J): 500 J')

#  ('Pulse (ps)', 10), ('Energy (J)', 900), ('Thick. (um)', varies))
temp(35431, None, 'blue', title='PJ-X 2; Energy (J): 900 J')
temp(35433, 35440, 'red', title='PJ-X 2; Energy (J): 900 J')
temp(35435, 35437, 'green', title='PJ-X 2; Energy (J): 900 J')


#  ('Pulse (ps)', 10), ('Energy (J)', 900), ('Thick. (um)', varies))
# labels = ['Energy (J)']
# temp(35427, None, 'blue', title='PJ-X 2; Pulse (ps): 10; Thick. (um): 0.5')
# temp(35429, None, 'red', title='PJ-X 2; Pulse (ps): 10; Thick. (um): 0.5')
# temp(35431, None, 'green', title='PJ-X 2; Pulse (ps): 10; Thick. (um): 0.5')


ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))



plt.show()