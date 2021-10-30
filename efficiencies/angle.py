import pickle
import numpy as np
from pathlib import Path
from lmfit.model import ModelResult
from matplotlib import pyplot as plt
from JSB_tools import fill_between


angel_x = []
angel_y = []
with open(Path(__file__).parent/'angel_3mm') as f:
    lines = f.readlines()
    for line in lines[1:]:
        x, _, y = list(map(float, line.split()))
        angel_x.append(x)
        angel_y.append(y)



with open(Path(__file__).parent/'3_17mm.pickle', 'rb') as f:
    model: ModelResult = pickle.load(f)


ergs = np.arange(60, 1500, 1)
y = model.eval(x=ergs)
yerr = model.eval_uncertainty(x=ergs)

ax = fill_between(ergs, y, yerr)
ax.plot(angel_x, angel_y)

plt.show()