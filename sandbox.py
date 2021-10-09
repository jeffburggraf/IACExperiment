# import plotly.graph_objects as go
# animals=['giraffes', 'orangutans', 'monkeys']
#
# fig = go.Figure(data=[
#     go.Bar(name='SF Zoo', x=animals, y=[20, 14, 23]),
#     go.Bar(name='LA Zoo', x=animals, y=[12, 18, 29])
# ])
# fig.update_layout(barmode='group')
#
# fig.show()
from pathlib import Path
import pandas
# from pandas import DataFrame
# import plotly.graph_objects as go
import numpy as np
a = np.array([1,2,3,4,5,6,])

p = Path(r"C:\Users\garag\OneDrive\Desktop\fakeSPE.Spe")
with open(p) as f:
    lines = f.readlines()
for l_i, line in enumerate(lines):
    if '@@@' in line:
        counts = np.zeros(16384, dtype=int)
        for i in range(0, 100, 10):
            counts[i] = 100
        line = '\n'.join(map(str, counts)) + '\n'
        lines[l_i] = line
print(p.with_name('fuckyou').with_suffix('.Spe'))
with open(p.with_name('fuckyou').with_suffix('.Spe'), 'w') as f:
    f.write(''.join(lines))

print(np.where(a>0))
