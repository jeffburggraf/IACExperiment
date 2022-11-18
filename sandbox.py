import matplotlib.pyplot as plt
import numpy as np


def get_mean_free_path(carrier_atoms, carrier_fracs, ff_lightQ, pressure=1):

    kb = 1.380649E-23  # J/K
    T = 270  # k
    ff = 'Kr' if ff_lightQ else 'Xe'

    ds = {'He': 260E-12, 'Ar': 340E-12, 'Xe': 396E-12, 'Kr': 360E-12}

    d_carrier = np.average([ds[a] for a in carrier_atoms], weights=carrier_fracs)

    d = np.mean([d_carrier, ds[ff]])

    pressure *= 1E5  # bar -> Pa

    out = kb * T / (np.sqrt(2) * np.pi * d ** 2 * pressure)
    return out


def get_D(atoms,carrier_fracs,  ff_lightQ, pressure=1):
    """

    Args:
        atoms:
        N: Number of nucleons of species which is diffusing
        pressure:

    Returns:

    """
    kb = 1.380649E-23  # J/K
    T = 270  # k

    N = 140 if ff_lightQ else 95

    lambda_ = get_mean_free_path(atoms, carrier_fracs, ff_lightQ, pressure=pressure)
    kg_per_u = 1.660539066E-27
    ke = 3/2 * kb * T

    mass = kg_per_u * N
    mean_vel = np.sqrt(2/mass * ke)
    return mean_vel * lambda_/3


print()

ts = np.linspace(0, 20, 100)

mean_sqr_displacments = {139: get_D(['He', 'Ar'], [1, 1], False),
                         95:  get_D(['He', 'Ar'], [1, 1], True)}


for k, d in mean_sqr_displacments.items():
    v = 100*np.sqrt(2 * ts * d)
    print(f"D = {d} for {k}")
    plt.plot(ts, v, label=k)
plt.legend()

plt.show()