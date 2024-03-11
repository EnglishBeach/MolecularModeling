from pathlib import Path
import os
import openmm
import openmm.app as openmm_app
import pandas as pd
from tools import base, simulator


def loop():
    for path in os.listdir('Find_D/systems'):
        sim_path = Path(f'Find_D/systems/{path}')
        substance, x, rho = path.split('_')
        result_dir = Path(f'Find_D/results/{sim_path.stem}')
        result_dir.mkdir(parents=True, exist_ok=True)

        print('#' * 30, f'{substance:10}: {int(x):<3}')
        for T in [15, 20, 25, 30, 40]:
            print('#' * 20, 'Tempetature: ', T)
            sim = simulator.Simulation(T=T)
            sim.load(sim_path)
            sim.equilibrate(20)
            sim.simulate(1000)
            df: pd.DataFrame = sim.get_data()
            df.to_csv(f'{result_dir}/{T}.csv')


loop()
