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
import pandas
import plotly.graph_objects as go
import numpy as np

# Create figure
fig = go.Figure()


xs = np.linspace(0, 10, 12)
ys = [i*xs**2 for i in range(len(xs))]
# print(ys)
# for k, v in {t: d for t, d in zip(xs, ys)}.items():
#     print(k, v)
df = pandas.DataFrame({t: d for t, d in zip(xs, ys)})
print(df)
df = df.melt()
print(df)
# Add traces, one for each slider step
for step in np.arange(0, 5, 0.1):
    fig.add_trace(
        go.Bar(
            visible=False,
            name="𝜈 = " + str(step),
            x=np.arange(0, 10, 0.01),
            y=np.sin(step * np.arange(0, 10, 0.01))))

# Make 10th trace visible
fig.data[10].visible = True

# Create and add slider
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
              {"title": "Slider switched to step: " + str(i)}],  # layout attribute
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)
# for s in steps:
#     print(s)

sliders = [dict(
    active=10,
    currentvalue={"prefix": "Frequency: "},
    pad={"t": 50},
    steps=steps
)]

fig.update_layout(
    sliders=sliders
)

# fig.show()