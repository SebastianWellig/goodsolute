import subprocess
import sys
from pathlib import Path

import pytest

ROOT_DIR = Path(__file__).resolve().parent.parent


@pytest.mark.parametrize(
    'case',
    [
        pytest.param(
            {
                'solute': 'tests/data/CH3OH.xyz',
                'solvent': 'tests/data/H2O.mol',
                'options': {'-p': 78.3, '-d': 0.997},
                'expected': {'dS': -0.0201528, 'dG': 6.00856},
            },
            id='methanol-in-water-RT',
        ),
        pytest.param(
            {
                'solute': 'tests/data/CH3OH.xyz',
                'solvent': 'tests/data/H2O.mol',
                'options': {'-p': 68.3, '-d': 0.988, '-t': 323.15},
                'expected': {'dS': -0.0202251, 'dG': 6.53575},
            },
            id='methanol-in-water-323.15K',
        ),
        pytest.param(
            {
                'solute': 'tests/data/Cl-.xyz',
                'solvent': 'tests/data/H2O.mol',
                'options': {'-p': 78.3, '-d': 0.997},
                'expected': {'dS': -0.0167807, 'dG': 5.00317},
            },
            id='chloride-in-water-RT',
        ),
    ],
)
def test_cli_examples(case):
    """Test examples with CLI."""
    command = [sys.executable, '-m', 'goodsolute', case['solute'], '-s', case['solvent']]
    for option, value in case['options'].items():
        command.extend([option, str(value)])
    result = subprocess.run(command, capture_output=True, text=True, check=True, cwd=ROOT_DIR)
    outputs = {}

    for line in result.stdout.splitlines():
        if line.startswith('Solvation correction to '):
            key, value = line.split('=', 1)
            outputs[key.strip()] = float(value.split(maxsplit=1)[0])

    expected = case['expected']
    assert outputs['Solvation correction to entropy'] == pytest.approx(expected['dS'], abs=5e-8)
    assert outputs['Solvation correction to Gibbs energy'] == pytest.approx(expected['dG'], abs=5e-6)
