from JSB_tools.PHITS_tools.PHITStoROOT import phits_to_root
from pathlib import Path
from JSB_tools import TBrowser, FileManager
import re
#
my_th_dict = {0: 0, 1: 2.5, 2: 5, 3: 10, 4: 20}

fman = FileManager(Path(__file__).parent, recreate=True)

for path in Path(__file__).parent.iterdir():
    if path.is_dir():
        if m := re.match('([0-9]+)_inp', path.name):
            my_th = my_th_dict[int(m.groups()[0])]
            ptrac_path = path/'PTRAC.txt'
            fman.add_path(ptrac_path.with_suffix('.root'), missing_ok=True, overwrite_ok=True, ptrac=True, my_thickness=my_th)
            phits_to_root(ptrac_path)

# TBrowser()
