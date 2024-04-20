import enum

from rdkit import Chem


class Mols(enum.Enum):
    butanol = 'BUT'
    ocm = 'OCM'
    dmag = 'DMA'
    peta = 'PET'


RD_MOLECULES = {
    Mols.ocm: Chem.MolFromSmiles("CC(C(OCCOC(OCCOCCOC(OCCOC(C(C)=C)=O)=O)=O)=O)=C"),
    Mols.dmag: Chem.MolFromSmiles("C=C(C)C(OCCOC(C(C)=C)=O)=O"),
    Mols.peta: Chem.MolFromSmiles("C=CC(OCC(COC(C=C)=O)(COC(C=C)=O)CO)=O"),
    Mols.butanol: Chem.MolFromSmiles("OCCCC"),
}

DATAS = {
    Mols.butanol: [
        (100, 700),
    ],
    Mols.dmag: [
        (0, 1069),
        # (11, 1069),
        # (22, 1050),
        # (35, 1029),
        # (39, 1019),
    ],
    Mols.ocm: [
        (0, 1520),
        # (12, 1580),
        # (22, 1569),
        # (31, 1550),
        # (38, 1539),
        # (44, 1530),
        # (49, 1510),
        # (54, 1489),
        # (58, 1490),
    ],
    Mols.peta: [
        (0, 1200),
        # (19, 1180),
        # (30, 1159),
        # (41, 1140),
        # (49, 1110),
    ],
}
