# goodsolute

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg?logo=python&logoColor=white)](https://github.com/terpy-tools/goodsolute/blob/HEAD/pyproject.toml)
[![Issues](https://img.shields.io/github/issues-raw/terpy-tools/goodsolute.svg?label=Issues)](https://github.com/terpy-tools/goodsolute/issues)
[![License](https://img.shields.io/badge/License-MIT-orange.svg)](https://github.com/terpy-tools/goodsolute/blob/HEAD/LICENSE)

CLI tool for computing [Garza's entropy corrections](#references) for solvation from geometry of solute and solvent molecules. Uses the scaled-particle theory $S_\epsilon$ for cavity corrections. Prints total entropy correction in kcal/(mol·K) and Gibbs energy correction in kcal/mol.

## Requirements

- Python 3.8 or newer.
- [RDKit](https://www.rdkit.org/)
- [cclib](https://cclib.github.io/)
- [SciPy](https://scipy.org/)

## Installation

Install directly from GitHub using `pip` and `git`:

```sh
python3 -m pip install git+https://github.com/terpy-tools/goodsolute.git
```

## Usage
### CLI

```sh
python3 -m goodsolute --help
```

| Argument | Description | Default |
| --- | --- | --- |
| `solute` | File with the coordinates for a solute molecule (any format supported by RDKit or cclib). | Required |
| `-s`<br />`--solvent` | File with the coordinates for a solvent molecule (any format supported by RDKit or cclib). | Required |
| `-t`<br />`--temperature` | Temperature [K]. | `298.15` |
| `-p`<br />`--permittivity` | Relative permittivity of the solvent. | Required |
| `-d`<br />`--density` | Density of the pure solvent [g/ml]. | Required |

Example usage:

```sh
# Compute corrections in water at the default temperature.
python3 -m goodsolute solute.log -s water.xyz -p 78.4 -d 0.997
```

Example output:
```log
--------------------------------------------------------------------------------
Garza's entropy corrections for solvation.

A. J. Garza, J. Chem. Theory Comput. 2019, 15, 3204-3214.
https://doi.org/10.1021/acs.jctc.9b00214

J. Ariai, U. Gellrich, Phys. Chem. Chem. Phys. 2023, 25, 14005-14015.
https://doi.org/10.1039/D3CP00970J

--------------------------------------------------------------------------------
Solvation correction to entropy      = -0.0152192 kcal/(mol*K)
Solvation correction to Gibbs energy = 4.53759 kcal/mol
--------------------------------------------------------------------------------
```

### As Python Library

Example usage:
```py
from rdkit import Chem
from rdkit.Chem import AllChem
from goodsolute.models import Solvent, Solute
from goodsolute.entropy_corrections import compute_total_correction

# Build solvent (H2O)
solvent_mol = Chem.MolFromSmiles('O')
solvent_mol = Chem.AddHs(solvent_mol)
AllChem.EmbedMolecule(solvent_mol, AllChem.ETKDG())

# Build solute (CH3OH)
solute_mol = Chem.MolFromSmiles('CO')
solute_mol = Chem.AddHs(solute_mol)
AllChem.EmbedMolecule(solute_mol, AllChem.ETKDG())

# Calculate correction for RDKit structures
solvent = Solvent.from_rdmol(solvent_mol, density=0.997, permittivity=78.3)
solute = Solute.from_rdmol(solute_mol)
correction = compute_total_correction(solvent, solute, temperature=298.15)
print(f'Entropy correction: {correction:6g} kcal/(mol*K)')
```

## Further Considerations
### Potential Double Counting from the Implicit Solvent Model

Some implicit solvation models, such as SMD, include non-electrostatic contributions $G_\mathrm{nes}$, which are fitted to reproduce $\Delta_\mathrm{solv} G$. [Ariai and Gellrich](#references) noted that there is a danger of double counting when applying Garza's corrections and therefore proposed that removing $G_\mathrm{nes}$ can be more accurate.

In the case of SMD, $G_\mathrm{nes}$ is usually reported as `SMD-CDS`, for instance in Gaussian log file:
```log
 SCF Done:  E(UMN15) =  -2160.62761811     A.U. after   20 cycles
            NFock= 20  Conv=0.90D-08     -V/T= 2.0041
 <Sx>= 0.0000 <Sy>= 0.0000 <Sz>= 0.5000 <S**2>= 0.7596 S= 0.5048
 <L.S>=  0.00000000000    
 SMD-CDS (non-electrostatic) energy       (kcal/mol) =      -2.99
 (included in total energy above)
 Annihilation of the first spin contaminant:
```

### Additional Corrections

The solvation correction $\Delta_\mathrm{solv} S$ only include terms for translational $\Delta_\mathrm{solv} S^\mathrm{trans}$, rotational $\Delta_\mathrm{solv} S^\mathrm{rot}$ and cavity entropy $S^\mathrm{cav}$. Additional corrections, such as quasi-harmonic and standard state corrections, have to be done separately, for instance, with [`goodvibes`](https://github.com/patonlab/goodvibes).



## References

1. [A. J. Garza, *J. Chem. Theory Comput.* **2019**, *15*, 3204-3214.](https://doi.org/10.1021/acs.jctc.9b00214)
2. [J. Ariai, U. Gellrich, *Phys. Chem. Chem. Phys.* **2023**, *25*, 14005-14015.](https://doi.org/10.1039/D3CP00970J)

## License

[MIT License](LICENSE).
