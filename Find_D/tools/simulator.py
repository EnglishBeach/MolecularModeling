import warnings
from pathlib import Path

import numpy as np
import openmm
import pandas as pd
from tqdm import tqdm as _tqdm

from .base import MolNames
from .boxer import Box

warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")


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
    check_freq=1000,
):

    # Integration options
    fs = openmm.unit.femtoseconds
    dt = dt * fs
    temperature = (T + 273) * openmm.unit.kelvin
    friction = 1 / openmm.unit.picosecond

    integrator = openmm.LangevinIntegrator(temperature, friction, dt)
    simulation = box.box_parametrized.to_openmm_simulation(integrator=integrator)

    equilibration = _tqdm(iterable=range(5))
    equilibration.set_description_str('Equilibration')
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.context.reinitialize(preserveState=True)
    for i in equilibration:
        simulation.step(1000)

    box.box_parametrized.to_pdb(
        boxes_path / Path(f"box_{box.substance.name}_{box.solvent_n}_{box.rho}.pdb"),
    )

    msdReporter = MSDReporter(check_freq, simulation, dt)
    simulation.reporters.append(msdReporter)
    simulation.currentStep = 0
    return simulation, msdReporter


def solve_D(T, box, result_dir: Path):
    simulation, reporter = create_simulation(
        box=box,
        T=T,
        boxes_path=result_dir,
    )
    product_cycle = _tqdm(iterable=range(500))
    product_cycle.set_description_str('Product      ')
    for i in product_cycle:
        simulation.step(1000)

    df: pd.DataFrame = reporter.df
    df.to_csv(f'{result_dir}/{T}.csv')
