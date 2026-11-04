"""Garza's entropy corrections for solvation.

A. J. Garza, J. Chem. Theory Comput. 2019, 15, 3204-3214.
https://doi.org/10.1021/acs.jctc.9b00214

J. Ariai, U. Gellrich, Phys. Chem. Chem. Phys. 2023, 25, 14005-14015.
https://doi.org/10.1039/D3CP00970J
"""

from math import log, pi

from scipy.constants import Boltzmann as k_B
from scipy.constants import R, angstrom, bar, calorie, kilo

from .objects import Solute, Solvent, SphericalCavity

STANDARD_PRESSURE = 1.0 * bar
R_KCAL = R / (kilo * calorie)


def compute_translational_correction(
    solvent: Solvent,
    solute: Solute,
    temperature: float,
) -> float:
    """Return correction to translational entropy [kcal/(mol*K)].

    :param solvent: Solvent instance.
    :param solute: Solute instance.
    :param temperature: Temperature [K].
    """
    cavity = SphericalCavity.from_molecules(solvent, solute)
    hopping_probability = max(solvent.free_volume ** (2 / 3) - solute.vdW_volume ** (2 / 3), 0.0) / (
        solvent.free_volume ** (2 / 3) + solvent.vdW_volume ** (2 / 3)
    )
    hopping_sites_count = 4 * cavity.volume ** (2 / 3) / (solvent.free_volume ** (2 / 3) + solvent.vdW_volume ** (2 / 3))
    accessible_cavities_count = 1 + hopping_sites_count * hopping_probability / (1 - hopping_probability)
    effective_volume = accessible_cavities_count * cavity.volume * angstrom**3
    reference_volume = k_B * temperature / STANDARD_PRESSURE
    return R_KCAL * log(effective_volume / reference_volume)


def compute_rotational_correction(
    solvent: Solvent,
    solute: Solute,
) -> float:
    """Return correction to rotational entropy [kcal/(mol*K)].

    :param solvent: Solvent instance.
    :param solute: Solute instance.
    """
    cavity = SphericalCavity.from_molecules(solvent, solute)
    min_radius_diff = (3 / 4 / pi) ** (2 / 3) * solvent.free_volume ** (1 / 3)
    radius_diff = max(cavity.radius - solute.radius_of_gyration, min_radius_diff)
    return 3 * R_KCAL * log(radius_diff / cavity.radius)


def compute_cavity_correction(
    solvent: Solvent,
    solute: Solute,
) -> float:
    """Return correction for cavity entropy based on S_epsilon approach [kcal/(mol*K)].

    :param solvent: Solvent instance.
    :param solute: Solute instance.
    """
    radius_ratio = (solute.vdW_volume / solvent.vdW_volume) ** (1 / 3)
    y = 3 / 4 / pi * (solvent.permittivity - 1) / (solvent.permittivity + 2)
    y_ratio = y / (1 - y)
    y_prime = 3 * y_ratio + 4.5 * y_ratio**2
    return -R_KCAL * (-log(1 - y) + 3 * y_ratio * radius_ratio + y_prime * radius_ratio**2)


def compute_total_coorection(
    solvent: Solvent,
    solute: Solute,
    temperature: float,
) -> float:
    """Return total correction to entropy [kcal/(mol*K)].

    :param solvent: Solvent instance.
    :param solute: Solute instance.
    :param temperature: Temperature [K].
    """
    return (
        compute_translational_correction(solvent, solute, temperature=temperature)
        + compute_rotational_correction(solvent, solute)
        + compute_cavity_correction(solvent, solute)
    )
