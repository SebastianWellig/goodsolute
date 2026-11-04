"""Classes for molecules and cavity properties."""

from math import pi

from rdkit import Chem
from rdkit.Chem.AllChem import ComputeMolVolume
from rdkit.Chem.Descriptors import MolWt
from rdkit.Chem.rdMolDescriptors import CalcRadiusOfGyration
from scipy.constants import Avogadro as N_A
from scipy.constants import angstrom, liter, milli


class Molecule:
    """Molecular properties.

    :param molar_mass: Molar mass [g/mol].
    :param vdW_volume: Van der Waals volume [A^3].
    """

    def __init__(self, molar_mass: float, vdW_volume: float):
        self.molar_mass = molar_mass
        self.vdW_volume = vdW_volume

    @classmethod
    def from_rdmol(cls, mol: Chem.Mol, **kwargs) -> 'Molecule':
        """Create instance from an RDKit molecule"""
        return cls(
            molar_mass=MolWt(mol),
            vdW_volume=ComputeMolVolume(mol),
            **kwargs,
        )


class Solvent(Molecule):
    """Solvent properties.

    :param molar_mass: Molar mass [g/mol].
    :param vdW_volume: Van der Waals volume [A^3].
    :param density: Density [g/mL].
    :param permittivity: Relative permittivity.
    """

    def __init__(
        self,
        molar_mass: float,
        vdW_volume: float,
        density: float,
        permittivity: float,
    ):
        super().__init__(molar_mass, vdW_volume)
        self.density = density
        self.permittivity = permittivity
        density_per_angstrom3 = density * angstrom**3 / (milli * liter)
        self.free_volume = molar_mass / N_A / density_per_angstrom3 - vdW_volume

    @classmethod
    def from_rdmol(cls, mol: Chem.Mol, density: float, permittivity: float) -> 'Solvent':
        """Create instance from an RDKit molecule"""
        return super().from_rdmol(mol, density=density, permittivity=permittivity)


class Solute(Molecule):
    """Solute properties.

    :param molar_mass: Molar mass [g/mol].
    :param vdW_volume: Van der Waals volume [A^3].
    :param radius_of_gyration: Radius of gyration [A].
    """

    def __init__(
        self,
        molar_mass: float,
        vdW_volume: float,
        radius_of_gyration: float,
    ):
        super().__init__(molar_mass, vdW_volume)
        self.radius_of_gyration = radius_of_gyration

    @classmethod
    def from_rdmol(cls, mol: Chem.Mol) -> 'Solute':
        """Create instance from an RDKit molecule"""
        return super().from_rdmol(mol, radius_of_gyration=CalcRadiusOfGyration(mol, useAtomicMasses=False))


class SphericalCavity:
    """Solvation cavity with derived radius [A].

    :param volume: Cavity volume [A^3].
    """

    def __init__(self, volume: float):
        self.volume = volume
        self.radius = (3 * self.volume / 4 / pi) ** (1 / 3)

    @classmethod
    def from_molecules(cls, solvent: Solvent, solute: Solute) -> 'SphericalCavity':
        """Create instance from solvent and solute instances"""
        volume = (solute.vdW_volume ** (1 / 3) + solvent.free_volume ** (1 / 3)) ** 3
        return cls(volume)
