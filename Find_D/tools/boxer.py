from pathlib import Path

from openff import interchange, toolkit, units
from openff.interchange.components import _packmol as packmol

from .base import RD_MOLECULES, MolNames

OFF_MOLECULES = {}
for molecule_type, rdkit_mol in RD_MOLECULES.items():
    mol = toolkit.Molecule.from_rdkit(rdkit_mol)
    OFF_MOLECULES[molecule_type] = mol
    mol.generate_conformers(n_conformers=1)
    mol.name = molecule_type.value

    for atom in mol.atoms:
        atom.metadata["residue_name"] = molecule_type.value
    mol.add_hierarchy_scheme(
        iterator_name="residue",
        uniqueness_criteria=["residue_name"],
    )


class Box:
    box = None
    ff = toolkit.ForceField("openff_unconstrained-2.1.0.offxml")

    def __init__(self, x, rho, substance: MolNames):
        self.solvent_n: int = int(x)
        self.substance_n: int = int((100 - x))
        self.substance: MolNames = substance
        self.rho = rho

    def __repr__(self) -> str:
        return f"<Box: {self.substance.name}= {self.substance_n} ({self.solvent_n}), rho={self.rho} mg/cm3>"

    def pack(self, tol=0.5):
        solvent = OFF_MOLECULES[MolNames.butanol]
        substance = OFF_MOLECULES[MolNames(self.substance)]
        if self.substance_n == 0:
            molecules = [solvent]
            n_molecules = [100]
        elif self.substance_n == 100:
            molecules = [substance]
            n_molecules = [100]
        else:
            molecules = [solvent, substance]
            n_molecules = [self.solvent_n, self.substance_n]
        kg = units.unit.kilogram
        m = units.unit.meter
        A = units.unit.angstrom
        self.box = packmol.pack_box(
            molecules=molecules,
            number_of_copies=n_molecules,
            mass_density=self.rho * kg / m**3,
            tolerance=tol * A,
            box_shape=packmol.UNIT_CUBE,
        )

    def parametrize(self):
        self.box_parametrized = interchange.Interchange.from_smirnoff(
            force_field=self.ff,
            topology=self.box,
        )

    def minimaze(self):
        self.box_parametrized.minimize()

    def save(self, path: Path):
        with open(path, "w") as file:
            box_j = self.box.to_json()
            data = {
                "solvent_n": self.solvent_n,
                "substance": self.substance.value,
                "rho": self.rho,
            }
            box_j = f"{data} ###" + box_j
            file.write(box_j)

    @classmethod
    def load(cls, path: Path):
        with open(path, "r") as file:
            box_j = file.read()
            data, box_j = box_j.split("###")
            data = eval(data)
            box = Box(
                x=int(data["solvent_n"]),
                rho=data["rho"],
                substance=MolNames(data["substance"][:3].upper()),
            )
            box.box = toolkit.Topology.from_json(box_j)

        return box

    def create_system(self):
        self.box_parametrized.to_openmm_topology()
        return self.box_parametrized.to_openmm_system()
