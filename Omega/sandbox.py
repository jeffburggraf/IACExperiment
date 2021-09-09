import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from JSB_tools import convolve_gauss2d




a = np.zeros((1000,1000))
a[499,499] = 10
o = convolve_gauss2d(a, 10)
print(np.sum(o))
plt.imshow(o)
plt.colorbar()
plt.show()
