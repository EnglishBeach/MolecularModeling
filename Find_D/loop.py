import os
from pathlib import Path

import openmm.app as openmm_app
import pandas as pd
import tools.simulator as simulator
from tqdm import tqdm as _tqdm


def loop():
    for path in os.listdir('./Find_D/systems'):
        top_path = Path(f'./Find_D/systems/{path}') / Path('top.pdb')
        system_path = Path(f'./Find_D/systems/{path}') / Path('system.xml')
        result_dir = Path(f'./Find_D/results/{path}')
        result_dir.mkdir(parents=True, exist_ok=True)

        substance, x, rho = path.split('_')
        print('#' * 30, f'{substance:10}: {int(x):<3}')

        for T in [15, 20, 25, 30, 40]:
