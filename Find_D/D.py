import simulator
from tqdm import tqdm
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import os


def get_T_data(T, box, result_dir: Path):
    simulation, reporter = simulator.create_simulation(
        box=box,
        T=T,
        boxes_path=result_dir,
    )
    product_cycle = tqdm(iterable=range(500))
    product_cycle.set_description_str('Product      ')
    for i in product_cycle:
        simulation.step(1000)

    df: pd.DataFrame = reporter.df
    df.to_csv(f'{result_dir}/{T}.csv')
    D_data = {}
    for mol_name in list(df.columns)[1:]:
        linear_model = linregress(df['Time'], df[mol_name])
        slope = linear_model.slope
        D = slope / 6 * 24 * 60 * 60
        D_data[mol_name] = D

    return D_data


def loop():
    box_data = []
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
            T_data.append(get_T_data(T, box, result_dir))
            print()

        box_data.append(pd.DataFrame(T_data))
    box_df: pd.DataFrame = pd.compat(box_data)
    box_df.to_csv('./Find_D/D_data.csv')


loop()
