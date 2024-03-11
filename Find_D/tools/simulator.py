import warnings
from pathlib import Path

import numpy as np
import openmm
import openmm.app as openmm_app
import pandas as pd
from tqdm import tqdm as _tqdm
from typing_extensions import Unpack

from .base import MolNames
from .boxer import Box

warnings.filterwarnings("ignore")
warnings.simplefilter("ignore")


class MSDReporter:

    def __init__(self, interval: float, simulation, dt: float):
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
    def is_selected(atom, residue: str):
        selected_atoms = [openmm.app.Element.getBySymbol(symbol) for symbol in ['O', 'C', 'N']]

        select_residue = atom.residue.name == residue
        select_atom = atom.element in selected_atoms
        return select_residue and select_atom

    def describeNextReport(self, simulation):
        steps = self.interval - simulation.currentStep % self.interval
        return (steps, True, False, False, False, False)

    def report(self, simulation, state):
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


class Simulation:

    def __init__(self, T, dt: float = 1):
        fs = openmm.unit.femtoseconds
        self.dt = dt * fs
        self.T = (T + 273) * openmm.unit.kelvin
        self.check_freq: int = 1000
        self.friction = 1 / openmm.unit.picosecond

    def load(self, work_path: str):
        work_path = Path(work_path)
        system_path = work_path / Path('system.xml')
        topology_path = work_path / Path('top.pdb')
        pdb = openmm_app.PDBFile(topology_path.as_posix())
        topology, positions = pdb.topology, pdb.positions

        with open(system_path.as_posix(), 'r') as system_file:
            system: openmm.System = openmm.XmlSerializer.deserialize(system_file.read())

        self.create(system=system, topology=topology, positions=positions)

    def create(self, system, topology, positions):
        # ADD forces
        system.addForce(openmm.CMMotionRemover())

        integrator = openmm.LangevinIntegrator(self.T, self.friction, self.dt)

        platform = openmm.Platform.getPlatformByName('CUDA')
        platformProperties = {'Precision': 'single'}
        platformProperties["DeviceIndex"] = "0"

        self.simulation = openmm_app.Simulation(
            topology=topology,
            system=system,
            integrator=integrator,
            platform=platform,
            platformProperties=platformProperties,
        )
        self.simulation.context.setPositions(positions)

    def equilibrate(self, steps: int = 1):
        equilibration_stage = _tqdm(iterable=range(steps))
        equilibration_stage.set_description_str('Equilibration')
        self.simulation.minimizeEnergy()
        self.simulation.context.setVelocitiesToTemperature(self.T)
        self.simulation.context.reinitialize(preserveState=True)
        for i in equilibration_stage:
            self.simulation.step(self.check_freq)

        self.simulation.currentStep = 0
        self.simulation.context.setVelocitiesToTemperature(self.T)
        self.msd_reporter = MSDReporter(self.check_freq, self.simulation, self.dt)
        self.simulation.reporters.append(self.msd_reporter)

    def simulate(self, cycles=10):
        product_cycle = _tqdm(iterable=range(cycles))
        product_cycle.set_description_str('Product      ')
        for i in product_cycle:
            self.simulation.step(1000)

    def get_data(self):
        return pd.DataFrame(self.msd_reporter.data)
