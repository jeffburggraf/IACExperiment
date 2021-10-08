import struct
from struct import unpack, calcsize
from mpant_reader import MPA
import matplotlib.pyplot as plt
from JSB_tools import mpl_hist
import numpy as np

def read_32(f):
    s = calcsize('I')
    _bytes = f.read(s)
    out = unpack('I', _bytes)
    out = f"{out:032b}"


with open('fuck') as f:
   print(sum(map(int, f.readlines())))

with open('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/friday/MCA/shot120.lst', 'rb') as f:
    s = 0
    pos = None
    while True:
        try:
            f.read(1).decode('ascii')
            pos = f.tell()
        except UnicodeDecodeError:
            break
    f.seek(pos)
    l = {}
    h = {}
    adc2_data = []

    waiting4adc2 = False

    for i in range(3000000
                   ):
        b = f.read(2)[::-1]

        # print(unpack('I', b))
        try:
            bin_ = f"{unpack('H', b)[0]:016b}"[::-1]
        except struct.error as e:
            print(e)
            break
        low_word = bin_[:16]
        high_word = bin_[16:]
        # print(low_word)
        # if waiting4adc2:
            # adc2_data.append(int(b[:16]))
            # waiting4adc2 = False

        if high_word == '1010111100000101':
            waiting4adc2 = True

        if low_word not in ['1111111111111111', '0000100000000000']:
            try:
                l[low_word] += 1
            except KeyError:
                l[low_word] = 1
            # l(int(low_word, 2))
        if high_word not in ['1111111111111111', '0000000010000000']:
            try:
                h[high_word] += 1
            except KeyError:
                h[high_word] = 1
                # print(low_word, high_word)
            # if :
        # print(low_word, high_word)
        if i <100:
            print(bin_, b.hex())
            # print(bin_[:len(bin_)//2], bin_[len(bin_)//2:], b[:len(b)//2].hex(), b[len(b)//2:].hex())
        # print()
    # print(f.read(2000).decode('ascii'))
    l = {k: l[k] for k in sorted(l.keys())}
    h = {k: h[k] for k in sorted(h.keys())}
    print(h)
    # plt.bar(range(len(l)), list(l.values()), align='center')
    # plt.xticks(range(len(l)), list(l.keys()), rotation='vertical')
    # plt.subplots_adjust(bottom=0.4)
    # plt.plot(l.values())
    # plt.plot(h.values())
    # plt.plot(l)
    values, bins = np.histogram(adc2_data, bins=np.arange(0, max(adc2_data)))
    mpl_hist(bins, values, label='adc2_data')
plt.show()
