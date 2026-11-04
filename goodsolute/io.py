"""Read coordinate files with RDKit and quantum-chemistry outputs with cclib.

Input geometries should include explicit hydrogen coordinates.
"""

from functools import partial
from pathlib import Path

from cclib.io import ccread
from cclib.parser.data import ccData
from rdkit import Chem

from .objects import Solute, Solvent


def _disable_implicit_atoms(mol: Chem.Mol):
    if mol is None:
        return None

    for atom in mol.GetAtoms():
        atom.SetNoImplicit(True)

    mol.UpdatePropertyCache(strict=False)
    return mol


RDREAD_MAP = {
    '.xyz': lambda p: _disable_implicit_atoms(Chem.MolFromXYZFile(p)),
    '.mol': partial(Chem.MolFromMolFile, removeHs=False),
    '.sdf': lambda p: Chem.SDMolSupplier(p, removeHs=False)[0],
    '.mol2': partial(Chem.MolFromMol2File, removeHs=False),
    '.pdb': partial(Chem.MolFromPDBFile, removeHs=False),
}


def ccData_to_rdmol(data: ccData, coords_idx=-1) -> Chem.Mol:
    """Converts cclib data to RDKit mol

    :param data: cclib data.
    :param coords_idx: Index of the coordinates to read from cclib data (last by default).
    """
    builder = Chem.RWMol()

    for atomic_num in data.atomnos:
        atom = Chem.Atom(int(atomic_num))
        builder.AddAtom(atom)

    mol = builder.GetMol()
    conformer = Chem.Conformer(len(data.atomnos))
    conformer.Set3D(True)

    for idx, pos in enumerate(data.atomcoords[coords_idx]):
        conformer.SetAtomPosition(idx, pos)

    mol.AddConformer(conformer)
    _disable_implicit_atoms(mol)
    return mol


def read_rdmol(mol_file: Path) -> Chem.Mol:
    """Read a molecule, using the final geometry of comp chem log files.

    :param mol_file: Path to chemical file format supported by RDKit or cclib.
    """
    mol_file = Path(mol_file)
    if not mol_file.is_file():
        raise FileNotFoundError(mol_file)

    rdread = RDREAD_MAP.get(mol_file.suffix.lower())

    if rdread is not None:
        mol = rdread(str(mol_file))
    else:
        data = ccread(str(mol_file))

        if data is None:
            raise ValueError(f'cclib could not read data from file: {mol_file}') from None

        mol = ccData_to_rdmol(data)

    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f'File has no atoms: {mol_file}')

    if mol.GetNumConformers() == 0:
        raise ValueError(f'File has no coordinates: {mol_file}')

    return mol


def read_solvent(
    mol_file: Path,
    density: float,
    permittivity: float,
) -> Solvent:
    return Solvent.from_rdmol(read_rdmol(mol_file), density=density, permittivity=permittivity)


def read_solute(mol_file: Path) -> Solute:
    return Solute.from_rdmol(read_rdmol(mol_file))
