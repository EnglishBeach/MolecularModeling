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



def post_loop():
    for box_file in os.listdir('./Find_D/boxes'):
        box_file = Path(box_file)
        probe_files = os.listdir('./Find_D/results')
        if box_file.stem in probe_files:
            continue

        result_dir = Path(f"./Find_D/results/{box_file.stem}")
        box = tools.boxer.Box.load(f"./Find_D/boxes/{box_file.name}")
        print(
            '#' * 20,
            f'{box.substance.name:10} but: {box.solvent_n:<3} | sub: {box.substance_n:<3}',
            '#' * 20,
        )
        print()
        print()
        result_dir.mkdir(parents=True, exist_ok=True)

        box.parametrize()
        box.minimaze()

        T_data = []
        for T in [15, 20, 25, 30, 40]:
            print('#' * 20, f'Temperature: {T:2} C')
            simulator.solve_D(T, box, result_dir)
            print()


post_loop()
