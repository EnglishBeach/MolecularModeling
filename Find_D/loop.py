import sys

sys.path.insert(0, 'Find_D')

import os
from pathlib import Path

import pandas as pd
import simulator


def loop():
    for box_file in os.listdir('./Find_D/boxes'):
        box_file = Path(box_file)
        result_dir = Path(f"./Find_D/results/{box_file.stem}")
        box = simulator.Box.load(f"./Find_D/boxes/{box_file.name}")
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
            print()
            simulator.solve_D(T, box, result_dir)
            print()


def post_loop():
    for box_file in os.listdir('./Find_D/boxes'):
        box_file = Path(box_file)
        probe_files = os.listdir('./Find_D/results')
        if box_file.stem in probe_files:
            continue

        result_dir = Path(f"./Find_D/results/{box_file.stem}")
        box = simulator.Box.load(f"./Find_D/boxes/{box_file.name}")
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
