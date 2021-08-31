from pathlib import Path
import re
from JSB_tools.list_reader import MaestroListFile


data_dir = Path(__file__).parent/'exp_data'


def get_maesto_list_shot_paths():
    out = {}
    for path in data_dir.iterdir():
        if path.is_dir() and re.match('.+day', path.name):
            for path in path.iterdir():
                if m := re.match(r'shot([0-9]+)\.Lis', path.name):
                    out[int(m.groups()[0])] = path
    out = {k:v for k, v in sorted(out.items(), key=lambda x: x[0])}

    return out


def get_mpant_mca_shot_paths():
    out = {}
    for path in data_dir.iterdir():
        if path.is_dir() and re.match('.+day', path.name):
            mca_path = path/'MCA'
            if not mca_path.exists():
                continue
            for path in mca_path.iterdir():
                if m := re.match(r'shot([0-9]+)\.mpa', path.name):
                    out[int(m.groups()[0])] = path
    out = {k:v for k, v in sorted(out.items(), key=lambda x: x[0])}
    return out

print(get_mpant_mca_shot_paths())
if __name__ == '__main__':
    for path in get_maesto_list_shot_paths().values():
        m = MaestroListFile(path)
        m.pickle()
