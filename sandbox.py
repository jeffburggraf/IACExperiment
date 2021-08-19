import numpy as np

a = np.arange(5)
def get_bin(x):

    return np.searchsorted(a, x, side="right") - 1

print(a)

print(get_bin(4))