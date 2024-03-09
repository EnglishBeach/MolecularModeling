import warnings
from pathlib import Path

import numpy as np
import openmm
import pandas as pd
from scipy.stats import linregress
from tqdm import tqdm

from .general import OFF_MOLECULES, RD_MOLECULES, Box, MolNames

warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")

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


class MSDReporter:

    def __init__(
        self,
        interval: float,
        simulation,
        dt: float,
    ):
        self.interval = interval
        self.dt = dt.value_in_unit(openmm.unit.second)
        self.ids = {}

        self.data = {'Time': np.array([])}
        self.atom_map = {}
        self.positions = {}
        self.start_positions = {}

        self.residues = list({residue.name for residue in simulation.topology.residues()})
        state0 = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
        atoms = list(simulation.topology.atoms())
        for residue_name in self.residues:
            atom_ids = [
                atoms.index(atom) for atom in atoms if MSDReporter.is_selected(atom, residue_name)
            ]

            self.atom_map[residue_name] = atom_ids

            positions = state0.getPositions(asNumpy=True).value_in_unit(openmm.unit.centimeter)
            self.start_positions[residue_name] = positions[atom_ids]

            self.data.update({residue_name: np.array([])})

    @staticmethod
    def is_selected(
        atom,
        residue: str,
    ):
        selected_atoms = [openmm.app.Element.getBySymbol(symbol) for symbol in ['O', 'C', 'N']]

        select_residue = atom.residue.name == residue
        select_atom = atom.element in selected_atoms
        return select_residue and select_atom

    def describeNextReport(
        self,
        simulation,
    ):
        steps = self.interval - simulation.currentStep % self.interval
        return (steps, True, False, False, False, False)

    def report(
        self,
        simulation,
        state,
    ):
        for residue in self.residues:
            atom_ids = self.atom_map[residue]
            positions = state.getPositions(asNumpy=True).value_in_unit(openmm.unit.centimeter)[
                atom_ids
            ]
            start_positions = self.start_positions[residue]

            msd = np.mean(np.linalg.norm((start_positions - positions), axis=1) ** 2)
            self.data[residue] = np.append(self.data[residue], msd)
        time = state.getStepCount() * self.dt
        self.data['Time'] = np.append(self.data['Time'], time)

    @property
    def df(self):
        return pd.DataFrame(self.data)


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


def get_D(T, box, result_dir: Path):
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
