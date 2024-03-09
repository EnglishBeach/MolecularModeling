import os
import warnings
from pathlib import Path

import openmm
import pandas as pd
import simulator
from scipy.stats import linregress
from tqdm import tqdm

from Find_D.simulator import Box, MSDReporter

warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")
from simulator import MolNames

DATAS = {
    MolNames.butanol: [
        (100, 700),
    ],
    MolNames.dmag: [
        (0, 1069),
        (11, 1069),
        (22, 1050),
        (35, 1029),
        (39, 1019),
    ],
    MolNames.ocm: [
        (0, 1720),
        (12, 1580),
        (22, 1569),
        (31, 1550),
        (38, 1539),
        (44, 1530),
        (49, 1510),
        (54, 1489),
        (58, 1490),
    ],
    MolNames.peta: [
        (0, 1200),
        (19, 1180),
        (30, 1159),
        (41, 1140),
        (49, 1110),
    ],
}


def create_simulation(
    box: Box,
    dt=1,
    T=25,
    boxes_path: Path = Path('.'),
    check_freq=100,
):

    # Integration options
    dt = dt * openmm.unit.femtoseconds  # simulation timestep
    temperature = (T + 273) * openmm.unit.kelvin  # simulation temperature
    friction = 1 / openmm.unit.picosecond  # friction constant

    integrator = openmm.LangevinIntegrator(temperature, friction, dt)
    simulation = box.box_parametrized.to_openmm_simulation(integrator=integrator)

    equilibration = tqdm(iterable=range(50))
    equilibration.set_description_str('Equilibration')
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.context.reinitialize(preserveState=True)
    for i in equilibration:
        simulation.step(100)

    box.box_parametrized.to_pdb(
        f"{boxes_path}/box_{box.substance.name}_{box.solvent_n}_{box.rho}.pdb",
    )

    msdReporter = MSDReporter(check_freq, simulation, dt)
    simulation.reporters.append(msdReporter)
    simulation.currentStep = 0
    return simulation, msdReporter


def get_T_data(T, box, result_dir: Path):
    simulation, reporter = create_simulation(
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
