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

with open('/Users/burggraf1/PycharmProjects/IACExperiment/exp_data/tuesday/MCA/shot2.lst', 'rb') as f:
    s = 0
    pos = None
    while True:
        try:
            f.read(1).decode('ascii')
            pos = f.tell()
        except UnicodeDecodeError:
            break
    f.seek(pos)

    def read_word():
        byte1 = f.read(1)
        byte2 = f.read(1)
        # convert the byte to an integer representation
        byte1 = ord(byte1)
        byte2 = ord(byte2)
        # now convert to string of 1s and 0s
        byte1 = bin(byte1)[2:].rjust(8, '0')
        byte2 = bin(byte2)[2:].rjust(8, '0')
        # now byte contains a string with 0s and 1s
        return byte2, byte1


    def hextobin(h):
        return bin(int(h, 16))[2:].zfill(len(h) * 4)

    for i in range(3000000):

        # if i < 100:
        #     print([read_word() for i in range(2)])
        if i < 100:
            outh = ""
            outb = ""
            for i in range(2):
                b1 = f.read(1)[::-1]
                b2 = f.read(1)[::-1]
                # b1.reverse()
                # b2.reverse()
                outh += f" {b2.hex()} {b1.hex()}"
                outb += f" {hextobin(b2.hex())} {hextobin(b1.hex())}"
                # print(b2.hex(), b1.hex(), hextobin(b2.hex()), hextobin(b1.hex()))
                # b = f.read(1)[::-1]
            print(outh)
            print(outb)
            print()
            # b1 = f.read(1)[::-1]
            # b2 = f.read(1)[::-1]
            # b1.reverse()
            # b2.reverse()
            # print(b2.hex(), b1.hex(), hextobin(b2.hex()), hextobin(b1.hex()))
        #
        # # print(unpack('I', b))
        # try:
        #     bin_ = f"{unpack('I', b)[0]:032b}"
        # except struct.error as e:
        #     print(e)
        #     break
        # low_word = bin_[16:]
        # high_word = bin_[:16]
        # if i < 100:
        #     print(low_word, high_word, b.hex())
plt.show()
