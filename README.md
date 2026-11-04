# goodsolute

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg?logo=python&logoColor=white)](https://github.com/terpy-tools/goodsolute/blob/HEAD/pyproject.toml)
[![Issues](https://img.shields.io/github/issues-raw/terpy-tools/goodsolute.svg?label=Issues)](https://github.com/terpy-tools/goodsolute/issues)
[![License](https://img.shields.io/badge/License-MIT-orange.svg)](https://github.com/terpy-tools/goodsolute/blob/HEAD/LICENSE)

CLI tool for computing [Garza's entropy corrections](#References) for solvation using the scaled-particle theory $S_\epsilon$ for cavity corrections from geometry of solute and solvent molecules. Prints total entropy correction in kcal/(mol·K) and Gibbs energy correction in kcal/mol.

The solvation correction $\Delta_\mathrm{solv} S$ only include terms for translational $\Delta_\mathrm{solv} S^\mathrm{trans}$, rotational $\Delta_\mathrm{solv} S^\mathrm{rot}$ and cavity entropy $S^\mathrm{cav}$. Additional corrections, such as quasi-harmonic and standard state corrections, have to be done separately, for instance, with [`goodvibes`](https://github.com/patonlab/goodvibes).

## Requirements

- Python 3.8 or newer.
- [RDKit](https://www.rdkit.org/) for molecular properties and van der Waals volumes.
- [cclib](https://cclib.github.io/) for reading quantum-chemistry output files.
- [SciPy](https://scipy.org/) for physical constants and unit conversions.

## Installation

Install directly from GitHub using `pip` and `git`:

```sh
python3 -m pip install git+https://github.com/terpy-tools/goodsolute.git
```

## Usage

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


## References

- [A. J. Garza, *J. Chem. Theory Comput.* **2019**, *15*, 3204-3214.](https://doi.org/10.1021/acs.jctc.9b00214)
- [J. Ariai, U. Gellrich, *Phys. Chem. Chem. Phys.* **2023**, *25*, 14005-14015.](https://doi.org/10.1039/D3CP00970J)

## License

[MIT License](LICENSE).
