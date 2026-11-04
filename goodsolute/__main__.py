from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

from .entropy_corrections import __doc__ as entropy_corrections_doc
from .entropy_corrections import compute_total_correction
from .io import read_solute, read_solvent


def main():
    parser = ArgumentParser(
        description="CLI tool for computing Garza's entropy corrections for solvation.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        'solute',
        type=Path,
        help='File with the coordinates for a solute molecule (any format supported by RDKit or cclib).',
    )
    parser.add_argument(
        '-s',
        '--solvent',
        type=Path,
        help='File with the coordinates for a solvent molecule (any format supported by RDKit or cclib).',
    )
    parser.add_argument(
        '-t',
        '--temperature',
        type=float,
        default=298.15,
        help='Temperature [K].',
    )
    parser.add_argument(
        '-p',
        '--permittivity',
        type=float,
        required=True,
        help='Relative permittivity of the solvent.',
    )
    parser.add_argument(
        '-d',
        '--density',
        type=float,
        required=True,
        help='Density of the pure solvent [g/ml].',
    )
    args = parser.parse_args()
    solvent = read_solvent(args.solvent, density=args.density, permittivity=args.permittivity)
    solute = read_solute(args.solute)
    entropy_corrections = compute_total_correction(solvent, solute, temperature=args.temperature)
    print(80 * '-')
    print(entropy_corrections_doc)
    print(80 * '-')
    print(f'Solvation correction to entropy      = {entropy_corrections:6g} kcal/(mol*K)')
    print(f'Solvation correction to Gibbs energy = {-entropy_corrections * args.temperature:6g} kcal/mol')
    print(80 * '-')


if __name__ == '__main__':
    main()
