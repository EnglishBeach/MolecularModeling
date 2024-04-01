import os
from pathlib import Path

from openff import interchange, toolkit, units
from openff.interchange.components import _packmol as packmol
from openff.units import Quantity, unit

from base import MolNames

from . import configs

MOLECULES = {}
for mol_path in os.listdir('Compounds'):
    mol_path = Path(mol_path)
    mol = toolkit.Molecule.from_file(f'Compounds/{mol_path}')
    MOLECULES[MolNames[mol_path.stem]] = mol

FF = toolkit.ForceField("openff_unconstrained-2.1.0.offxml")


class Simulation:
    mpi = 4
    total_n = 200

    def __init__(self, workdir: Path, T, substance, x, rho):
        self.rho = rho
        self.x = x
        self.T = T + 273

        if x == 1.0:
            self.molecules = [MolNames.butanol]
            self.n_molecules = [self.total_n]
        elif x == 0.0:
            self.molecules = [MolNames[substance]]
            self.n_molecules = [self.total_n]
        else:
            self.molecules = [MolNames.butanol, MolNames[substance]]
            self.n_molecules = [
                int(self.x * self.total_n),
                self.total_n - int(self.x * self.total_n),
            ]

        self.WORK = workdir
        self.create_workspace()

        configs_path = workdir / '0configs'
        self.top = None
        self.RAW = workdir / '1box'

        self.EM = workdir / '2EM'
        self.em_config = configs.em(configs_path)

        self.NVT = workdir / '3NVT'
        self.nvt_config = configs.nvt(configs_path, T=self.T, compounds=self.molecules)

        self.NPT = workdir / '4NPT'
        self.npt_config = configs.npt(configs_path, T=self.T, compounds=self.molecules)

        self.MD = workdir / '5MD'
        self.md_config = configs.md(configs_path, T=self.T, compounds=self.molecules)

    def create_workspace(self):
        for stage in ['0configs', '1box', '2EM', '3NVT', '4NPT', '5MD']:
            folder = self.WORK / stage
            folder.mkdir(parents=True, exist_ok=True)

    def create_box(self):
        g = units.unit.gram
        cm = units.unit.centimeter
        A = units.unit.angstrom

        print('Pack')
        box = packmol.pack_box(
            molecules=[MOLECULES[mol] for mol in self.molecules],
            number_of_copies=self.n_molecules,
            mass_density=self.rho * g / cm**3,
            tolerance=0.1 * A,
            box_shape=packmol.UNIT_CUBE,
        )
        print('Parametrize')
        inter = interchange.Interchange.from_smirnoff(
            force_field=FF,
            topology=box,
            charge_from_molecules=[MOLECULES[mol] for mol in self.molecules],
        )
        print('Minimaze')
        inter.minimize(force_tolerance=Quantity(1.0, unit.kilojoule_per_mole / unit.nanometer))
        inter.to_pdb(self.RAW / 'box.pdb')
        inter.to_gromacs((self.RAW / 'box').as_posix())
        self.top = self.RAW / 'box.top'

    def center(self):
        print('*' * 20, ' Centering')
        command = f"""
gmx -quiet  editconf \
-f {self.RAW}/box.gro \
-o {self.RAW}/box.gro \
-c \
-bt cubic
"""
        res = os.system(command)
        if res:
            raise NotADirectoryError()

    def em(self):
        print('*' * 20, ' EM')
        command = f"""
gmx -quiet grompp \
-f  {self.em_config} \
-p {self.top} \
-c {self.RAW}/box.gro \
-maxwarn 5 \
-o {self.EM}/system.tpr \
-po {self.EM}/config.mdp

gmx -quiet mdrun \
-s {self.EM}/system.tpr \
-c {self.EM}/box.gro \
-deffnm {self.EM}/em \
-v \
-pin on \
-ntmpi {self.mpi}
"""
        res = os.system(command)
        if res:
            raise NotADirectoryError()

    def nvt(self):
        print('*' * 20, ' NVT')
        command = f"""
gmx -quiet grompp \
-f  {self.nvt_config} \
-p {self.top} \
-c {self.EM}/box.gro \
-maxwarn 5 \
-o {self.NVT}/system.tpr \
-po {self.NVT}/config.mdp

gmx -quiet mdrun \
-s {self.NVT}/system.tpr \
-c {self.NVT}/box.gro \
-deffnm {self.NVT}/nvt \
-v \
-pin on \
-ntmpi {self.mpi}
"""
        res = os.system(command)
        if res:
            raise NotADirectoryError()

    def npt(self):
        print('*' * 20, ' NPT')
        command = f"""
gmx -quiet grompp \
-f  {self.npt_config} \
-p {self.top} \
-c {self.NVT}/box.gro \
-maxwarn 5 \
-o {self.NPT}/system.tpr \
-po {self.NPT}/config.mdp

gmx -quiet mdrun \
-s {self.NPT}/system.tpr \
-c {self.NPT}/box.gro \
-deffnm {self.NPT}/npt \
-v \
-pin on \
-ntmpi {self.mpi}
    """
        res = os.system(command)
        if res:
            raise NotADirectoryError()

    def md(self):
        print('*' * 20, ' MD')
        command = f"""
gmx -quiet grompp \
-f  {self.md_config} \
-p {self.top} \
-c {self.NPT}/box.gro \
-maxwarn 5 \
-o {self.MD}/system.tpr \
-po {self.MD}/config.mdp

gmx -quiet mdrun \
-s {self.MD}/system.tpr \
-c {self.MD}/box.gro \
-deffnm {self.MD}/md \
-v \
-pin on \
-ntmpi {self.mpi}
    """
        res = os.system(command)
        if res:
            raise NotADirectoryError()

    def summarize(self):
        print('*' * 20, ' Resulting')
        command = f"""
gmx -quiet select \
-s {self.RAW}/box.gro \
-select 0 \
-on {self.MD}/index.ndx
"""
        res = os.system(command)
        if res:
            raise NotADirectoryError()

        command = f"""
gmx -quiet trjconv \
-f {self.MD}/md.xtc \
-s {self.MD}/system.tpr \
-n {self.MD}/index.ndx \
-center \
-pbc nojump \
-o {self.WORK}/traj.xtc
"""
        res = os.system(command)
        if res:
            raise NotADirectoryError()

        command = f"""
cp {self.MD}/box.gro {self.WORK}/box.gro
"""
        res = os.system(command)
        if res:
            raise NotADirectoryError()
