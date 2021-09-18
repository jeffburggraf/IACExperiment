from JSB_tools import isinstance_generic, mpl_hist
import numpy as np
from matplotlib import pyplot as plt
import time

#========================
n = 60000
hl = 200
bg_rel_n_points = 1
bg_scale = 1/3
tmax = 100
#  ===================================
sig_points = np.random.exponential(hl/np.log(2),  n)
sig, bins = np.histogram(sig_points, bins='auto')
bg_points = np.random.exponential(hl/2/np.log(2), int(n*bg_rel_n_points/bg_scale))
# bg_points = np.random.uniform(0, max(sig_points*1.5), int(n*bg_rel_n_points/bg_scale))
bg, _ = np.histogram(bg_points, bins=bins, weights=bg_scale*np.ones_like(bg_points))

tot_points = np.concatenate([sig_points, bg_points])
tot = sig+bg
ax = mpl_hist(bins, sig, label='sig')
mpl_hist(bins, tot, ax=ax, label='tot')
mpl_hist(bins, bg, ax=ax, label='bg')
# ax.legend()
tot_points = tot_points[np.where(tot_points < tmax)]
bg_points = bg_points[np.where(bg_points < tmax)]

lambda_est = (len(tot_points) - len(bg_points)) / (np.sum(tot_points) - np.sum(bg_points))
hl_est = np.log(2) / lambda_est
print('init: ', hl_est)

data_len = len(tot_points) - len(bg_points)
data_sum = np.sum(tot_points) - np.sum(bg_points)


i = 0


def iterate(prev_lambda=None):
    # global i
    # i += 1
    if prev_lambda is None:
        return iterate(data_len/data_sum)
    else:
        corr = -data_len * (tmax / (1 - np.e ** (tmax * prev_lambda)))
        # print(corr, data_len, data_sum)
        next_lambda = data_len/(data_sum + corr)
        # print(np.log(2)/next_lambda)
        if np.isclose(prev_lambda, next_lambda, rtol=0.000001):
            return next_lambda
        else:
            global cor
            return iterate(next_lambda)


t = time.time()
est = np.log(2)/iterate()
# print(corr, )
print(time.time() - t, 'n loops: ', i)
print(est)
plt.show()
