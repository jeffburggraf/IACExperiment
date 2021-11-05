from FFStopping import GetPoints
import numpy as np
from analysis import Shot


def get_fraction(shot: Shot, FF: str):
    if shot.foil_pos == 'upstream':
        def cut_up(x, y, z):
            if not np.linalg.norm([x, y]) <= 2:
                return False
            if not z < 1:
                return False
            return True

        def cut_down(x, y, z):
            if not np.linalg.norm([x, y]) <= 2:
                return False
            if not z < 5:
                return False
            return True

        if shot.ar_flow > 0 and shot.he_flow == 0 and shot.pressure > 1.2:
            points_debug = GetPoints(['Ar'], [1], 1.3)
            points_up = GetPoints(['Ar'], [1], 1.3, cut=cut_up)
            points_down = GetPoints(['Ar'], [1], 1.3, cut=cut_down)

        else:
            assert False
    else:
        assert False
    points_debug.plotz(FF)
    points_debug.plotr(FF)
    return points_up.percent_make_cut(FF), points_down.percent_make_cut(FF)