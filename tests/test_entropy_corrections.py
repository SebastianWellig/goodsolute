import pytest

from goodsolute.entropy_corrections import (
    compute_cavity_correction,
    compute_rotational_correction,
    compute_translational_correction,
)
from goodsolute.objects import Solute, Solvent, SphericalCavity


@pytest.mark.parametrize(
    'case',
    [
        pytest.param(
            {
                'solvent': Solvent(molar_mass=78.112, vdW_volume=78.857, density=0.8765, permittivity=2.2825),
                'solute': Solute(molar_mass=96.101, vdW_volume=82.82649, radius_of_gyration=2.052169),
                'trans_kwargs': {},
                'expected': {'V_cav': 606.164, 'dS_trans': -0.008382, 'dS_rot': -0.002955, 'S_cav': -0.001143},
            },
            id='Figure S4',
        ),
        pytest.param(
            {
                'solvent': Solvent(molar_mass=86.036, vdW_volume=80, density=0.968, permittivity=2),
                'solute': Solute(molar_mass=95.101, vdW_volume=82.82649, radius_of_gyration=2.052169),
                'trans_kwargs': {'temperature': 310},
                'expected': {'V_cav': 599.598, 'dS_trans': -0.008481, 'dS_rot': -0.002969, 'S_cav': -0.000929},
            },
            id='Figure S6',
        ),
    ],
)
def test_gellrich_numerical_examples(case):
    """Test numerical examples from SI of Phys. Chem. Chem. Phys. 2023, 25, 14005-14015."""
    solvent = case['solvent']
    solute = case['solute']
    expected = case['expected']
    cavity = SphericalCavity.from_molecules(solvent, solute)
    assert cavity.volume == pytest.approx(expected['V_cav'], abs=5e-4)
    assert compute_translational_correction(solvent, solute, **case['trans_kwargs']) == pytest.approx(expected['dS_trans'], abs=5e-7)
    assert compute_rotational_correction(solvent, solute) == pytest.approx(expected['dS_rot'], abs=5e-7)
    assert compute_cavity_correction(solvent, solute) == pytest.approx(expected['S_cav'], abs=5e-7)
