import warnings
from pathlib import Path

import numpy as np
import openmm
import openmm.app as openmm_app
import pandas as pd
from tqdm import tqdm as _tqdm
from openmm import unit

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
        state0 = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
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

    def __init__(
        self,
        T,
        dt_fs: float = 1,
        check_freq: int = 1000,
    ):
        fs = openmm.unit.femtoseconds
        ps = openmm.unit.picosecond
        K = openmm.unit.kelvin

        self.dt = dt_fs * fs
        self.T = (T + 273) * K
        self.check_freq = check_freq
        self.friction = 1 / ps
        self.reporter = {}

    def load(self, work_path: str):
        work_path = Path(work_path)
        system_path = work_path / 'system.xml'
        topology_path = work_path / 'top.pdb'
        pdb = openmm_app.PDBFile(topology_path.as_posix())
        topology, positions = pdb.topology, pdb.positions

        with open(system_path.as_posix(), 'r') as system_file:
            system: openmm.System = openmm.XmlSerializer.deserialize(system_file.read())

        self.create(system=system, topology=topology, positions=positions)

    @property
    def positions(self):
        return self.simulation.context.getState(
            getPositions=True, enforcePeriodicBox=True
        ).getPositions()

    def create(self, system, topology, positions):
        atmos = unit.atmospheres
        # ADD forces
        system.addForce(openmm.CMMotionRemover())
        system.addForce(openmm.MonteCarloBarostat(1 * atmos, self.T, 25))

        integrator = openmm.LangevinIntegrator(self.T, self.friction, self.dt)

        platform = openmm.Platform.getPlatformByName('CUDA')
        platformProperties = {}
        platformProperties['Precision'] = 'mixed'
        platformProperties['UseBlockingSync'] = 'true'

        self.simulation = openmm_app.Simulation(
            topology=topology,
            system=system,
            integrator=integrator,
            platform=platform,
            platformProperties=platformProperties,
        )
        self.simulation.context.setPositions(positions)

    def equilibrate(
        self,
        minimized_path: str,
        steps: int = 10,
        tolerance: float = 2,
    ):

        status = _tqdm(iterable=range(steps))
        status.set_description_str('Energy minimization')
        self.simulation.minimizeEnergy(
            maxIterations=self.check_freq * 10,
            tolerance=tolerance * unit.kilojoules_per_mole / unit.nanometer,
        )

        status.set_description_str('NVT/NPT minimisation')
        self.simulation.context.setVelocitiesToTemperature(self.T)
        self.simulation.context.reinitialize(preserveState=True)
        for i in status:
            self.simulation.step(self.check_freq)

        self.simulation.currentStep = 0
        self.simulation.context.setVelocitiesToTemperature(self.T)

        minimized_path: Path = Path(minimized_path)
        minimized_path.mkdir(exist_ok=True, parents=True)
        with open(minimized_path / 'Eq.pdb', 'w') as outfile:
            openmm_app.PDBFile.writeFile(
                self.simulation.topology,
                self.positions,
                file=outfile,
                keepIds=True,
            )

    def start_product(self, out_path: str):
        # ADD another forces
        out_path: Path = Path(out_path)
        out_path.mkdir(exist_ok=True, parents=True)

        self.reporter['msd'] = MSDReporter(
            self.check_freq,
            self.simulation,
            self.dt,
        )

        self.reporter['xtc'] = openmm_app.XTCReporter(
            file=(out_path / 'traj.xtc').as_posix(),
            reportInterval=self.check_freq,
            enforcePeriodicBox=True,
        )
        self.simulation.reporters.extend(list(self.reporter.values()))

    def simulate(self, cycles=10):
        status = _tqdm(iterable=range(cycles))
        status.set_description_str('Product')
        for i in status:
            self.simulation.step(self.check_freq)

    def get_msd(self):
        return pd.DataFrame(self.reporter['msd'].data)
