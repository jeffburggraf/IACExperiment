from __future__ import annotations
from datetime import datetime
import re
from JSB_tools import Nuclide
from typing import Union
from uncertainties import ufloat


class CalSource:
    instances = {}

    @classmethod
    def get_source(cls, serial) -> CalSource:
        try:
            return cls.instances[str(serial).lower()]
        except KeyError:
            assert False, f'No source entered for serial, "{serial}"\n. Options:\n{list(CalSource.instances.keys())}'

    def __init__(self, name, serial_num: Union[str, int], activity, ref_date: datetime, unit='uci', ):
        unit = unit.lower().strip()
        c = 1
        if 'ci' in unit:
            c *= 3.7E10
        if re.match("^u.{2}$", unit):
            c *= 1E-6
        elif re.match("^n.{2}$", unit):
            c *= 1E-9
        elif re.match("^m.{2}$", unit):
            c *= 1E-3
        elif re.match("^k.{2}$", unit):
            c *= 1E3

        self.name = name
        self.nuclide = Nuclide.from_symbol(name)
        self.serial_num = str(serial_num).lower()
        self._ref_activity = activity*c
        self.ref_date = ref_date

        CalSource.instances[self.serial_num] = self

    def get_activity(self, date=datetime.now()):
        dt = (date - self.ref_date).total_seconds()
        hl = self.nuclide.half_life
        correction = 0.5 ** (dt / hl)
        return self._ref_activity * correction

    def get_n_decays(self, duration, start_date=datetime.now()):
        dt = (start_date - self.ref_date).total_seconds()
        hl = self.nuclide.half_life
        n_nuclides_ref = self._ref_activity/self.nuclide.decay_rate
        n_nuclides_begin = n_nuclides_ref*0.5**(dt/hl)

        percent_decay = 1 - 0.5**(duration/hl)
        # print(n_nuclides_ref, percent_decay)
        return n_nuclides_begin*percent_decay

    def __repr__(self):
        return f'{self.name}; {self.get_activity():.3E} Bq'


CalSource("Y88", serial_num=190607000, activity=433, ref_date=datetime(2019, 7, 1), unit='kBq')
CalSource("Na22", serial_num=129742, activity=1.146, ref_date=datetime(2008, 7, 1), unit='uCi')
CalSource("Cs137", serial_num=129792, activity=92.11, ref_date=datetime(2008, 7, 1), unit='nCi')
CalSource("Cd109", serial_num=129757, activity=10.40, ref_date=datetime(2008, 7, 1), unit='uCi')
CalSource("Mn54", serial_num="J4-348", activity=9.882, ref_date=datetime(2012, 9, 1), unit='uCi')
CalSource("Co57", serial_num="K4-895", activity=1.135, ref_date=datetime(2013, 7, 15), unit='uCi')
