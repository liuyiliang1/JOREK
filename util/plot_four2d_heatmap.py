#!/usr/bin/env python3
"""
Read four2d output and plot heatmap (m vs psi_N) of Fourier amplitudes.
Usage: python plot_four2d_heatmap.py <four2d_output_dir> [options]

Output format from postproc four2d:
  - Files: {step_prefix}_absolute_value_n{nnn}.dat / _complex_phase_n{nnn}.dat
  - Header: # Psi_N  var1  var2 ...
  - Blocks: # {type} for m/n={m:+d}/{n:+d}
  - Data: psi_N, var1_amplitude, var2_amplitude, ...
"""

import numpy as np
import matplotlib.pyplot as plt
import re
import os
import glob
import argparse
from pathlib import Path


def find_four2d_file(data_dir, step, n_mode, output_type='absolute_value'):
    """Find the four2d output file matching given step, n, and output type."""
    patterns = [
        f"*_s{step:05d}_{output_type}_n{n_mode:03d}.dat",
        f"*s{step:05d}*{output_type}*n{n_mode:03d}*.dat",
    ]
    for pat in patterns:
        matches = glob.glob(os.path.join(data_dir, pat))
        if matches:
            return matches[0]
    return None


def parse_four2d_file(filepath):
    """Parse a single four2d output file, return dict of {m: array} and variable names list."""
    blocks = {}
    var_names = None
    current_m = None
    current_data = []

    with open(filepath, 'r') as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith('#') and 'm/n=' in line:
            if current_m is not None and current_data:
                blocks[current_m] = np.array(current_data)

            match = re.search(r'm/n\s*=\s*([+-]?\d+)\s*/\s*([+-]?\d+)', line)
            if match:
                current_m = int(match.group(1))
            current_data = []

        elif line.startswith('#') and 'Psi_N' in line and var_names is None:
            var_names = line.lstrip('#').split()

        elif not line.startswith('#'):
            vals = [float(x) for x in line.split()]
            current_data.append(vals)

    if current_m is not None and current_data:
        blocks[current_m] = np.array(current_data)

    return blocks, var_names


def list_available(data_dir):
    """Print available steps, n modes, variables, and output types in the directory."""
    files = sorted(glob.glob(os.path.join(data_dir, '*.dat')))
    if not files:
        print("No .dat files found.")
        return

    steps = set()
    n_modes = set()
    var_sets = {}
    out_types = set()

    for f in files:
        basename = os.path.basename(f)
        # Try to parse step, output_type, n from filename
        match = re.search(r's(\d{5})', basename)
        if match:
            steps.add(int(match.group(1)))

        match = re.search(r'n(\d{3})', basename)
        if match:
            n_modes.add(int(match.group(1)))

        for ot in ['absolute_value', 'complex_phase', 'real_part', 'imaginary_part']:
            if ot in basename:
                out_types.add(ot)

    # Read one file to get variable names
    try:
        _, var_names = parse_four2d_file(files[0])
        if var_names:
            var_names = var_names[1:]  # Skip Psi_N (coordinate)
            var_sets[os.path.basename(files[0])] = var_names
    except Exception:
        pass

    print(f"Found {len(files)} .dat files\n")
    print(f"Steps:     {sorted(steps)}")
    print(f"n modes:   {sorted(n_modes)}")
    print(f"Types:     {sorted(out_types)}")
    if var_names:
        print(f"Variables: {var_names}")
    else:
        print("Variables: (could not parse)")


def plot_heatmap(data_dir, step=1600, n_mode=0, var_name=None,
                 output_type='absolute_value', vmax=None, cmap='jet',
                 mrange=None, psirange=None):
    """Read four2d files and plot heatmap. var_name resolved automatically from file header."""

    filepath = find_four2d_file(data_dir, step, n_mode, output_type)
    if filepath is None:
        # List what's available
        print(f"No file found for step={step}, n={n_mode}, type={output_type}\n")
        list_available(data_dir)
        raise FileNotFoundError("Use --list to see all available files.")

    print(f"Reading: {filepath}")

    blocks, var_names = parse_four2d_file(filepath)

    if var_names is None:
        raise ValueError("Could not parse variable names from file header.")

    # Build dynamic var_map from header (skip column 0 = Psi_N)
    var_map = {name: idx for idx, name in enumerate(var_names)}

    # Resolve variable
    if var_name is None:
        non_coord = [v for v in var_names if v != 'Psi_N']
        var_name = non_coord[0] if non_coord else var_names[0]
        print(f"No --var specified, using first variable: {var_name}")

    if var_name not in var_map:
        raise ValueError(
            f"Variable '{var_name}' not found. Available: {list(var_map.keys())}"
        )

    var_index = var_map[var_name]

    print(f"Variables: {var_names}")
    print(f"Selected:  {var_name} (index {var_index})")
    print(f"Blocks found: {len(blocks)} (m from {min(blocks.keys())} to {max(blocks.keys())})")

    # Build the data matrix: rows = psi_N, columns = m
    m_sorted = sorted(blocks.keys())

    first_block = blocks[m_sorted[0]]
    psi_N = first_block[:, 0]
    n_psi = len(psi_N)
    n_m = len(m_sorted)

    data_matrix = np.zeros((n_psi, n_m))
    for j, m_val in enumerate(m_sorted):
        blk = blocks[m_val]
        data_matrix[:, j] = np.abs(blk[:, var_index])

    # ---- Crop m range ----
    if mrange is not None:
        m_min, m_max = mrange
        m_mask = [(m_val >= m_min) and (m_val <= m_max) for m_val in m_sorted]
        m_sorted = [m for i, m in enumerate(m_sorted) if m_mask[i]]
        data_matrix = data_matrix[:, m_mask]
        print(f"Cropped m to [{m_min}, {m_max}]: {len(m_sorted)} modes")

    # ---- Crop psi_N range ----
    if psirange is not None:
        psi_min, psi_max = psirange
        psi_mask = (psi_N >= psi_min) & (psi_N <= psi_max)
        psi_N = psi_N[psi_mask]
        data_matrix = data_matrix[psi_mask, :]
        print(f"Cropped psi_N to [{psi_min}, {psi_max}]: {len(psi_N)} points")

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(12, 6))

    extent = [m_sorted[0], m_sorted[-1], psi_N[0], psi_N[-1]]

    im = ax.imshow(data_matrix, aspect='auto', origin='lower',
                   extent=extent, cmap=plt.get_cmap(cmap),
                   vmin=0, vmax=vmax, interpolation='bicubic')

    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(var_name, fontsize=14)

    ax.set_xlabel('m (poloidal mode)', fontsize=14)
    ax.set_ylabel(r'$\psi_N$', fontsize=14)
    ax.set_title(f'{var_name}  |  n = {n_mode}  |  step = {step}  |  {output_type}',
                 fontsize=14)

    ax.tick_params(labelsize=12)
    plt.tight_layout()

    outname = f"{var_name}_n{n_mode:03d}_s{step:05d}_{output_type}.png"
    outpath = os.path.join(data_dir, outname)
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"Saved: {outpath}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Plot four2d Fourier heatmap (m vs psi_N). '
                    'Variable names are auto-detected from file header.')
    parser.add_argument('data_dir', help='Directory containing four2d .dat files')
    parser.add_argument('--step', type=int, default=1600, help='Time step')
    parser.add_argument('--n', type=int, default=0, help='Toroidal mode n')
    parser.add_argument('--var', type=str, default=None,
                        help='Variable name (auto-detected if omitted)')
    parser.add_argument('--type', type=str, default='absolute_value',
                        choices=['absolute_value', 'complex_phase'],
                        help='Output type')
    parser.add_argument('--vmax', type=float, default=None,
                        help='Max value for colorbar (auto if omitted)')
    parser.add_argument('--cmap', type=str, default='jet',
                        help='Colormap')
    parser.add_argument('--mmin', type=int, default=None,
                        help='Min poloidal mode m for cropping')
    parser.add_argument('--mmax', type=int, default=None,
                        help='Max poloidal mode m for cropping')
    parser.add_argument('--psimin', type=float, default=None,
                        help='Min psi_N for cropping')
    parser.add_argument('--psimax', type=float, default=None,
                        help='Max psi_N for cropping')
    parser.add_argument('-l', '--list', action='store_true',
                        help='List available steps, n modes, variables, and exit')

    args = parser.parse_args()

    if args.list:
        list_available(args.data_dir)
        return

    mrange = None
    if args.mmin is not None or args.mmax is not None:
        mrange = (args.mmin if args.mmin is not None else -999,
                  args.mmax if args.mmax is not None else 999)

    psirange = None
    if args.psimin is not None or args.psimax is not None:
        psirange = (args.psimin if args.psimin is not None else 0.0,
                    args.psimax if args.psimax is not None else 1.0)

    plot_heatmap(
        data_dir=args.data_dir,
        step=args.step,
        n_mode=args.n,
        var_name=args.var,
        output_type=args.type,
        vmax=args.vmax,
        cmap=args.cmap,
        mrange=mrange,
        psirange=psirange
    )


if __name__ == '__main__':
    main()
