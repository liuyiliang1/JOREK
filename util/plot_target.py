#!/usr/bin/env python3
"""
Plot ASCII outputs of the jorek2_target2vtkn diagnostic.

Two file types are auto-detected from the header:

1. target_surface2d.dat       -- 2D surface map (midplane2d-style layout):
   one '# Column m of n_plane  phi(deg)=...' block per toroidal angle,
   n_points rows per column ordered along the boundary, each row holding
   all values (es23.15).
   Columns: s  R  Z  psi  Psi_N  rho  density  Ti  Te  Vpar  angle  nv
            nvT_gam  KparT  bnd
   Default: colormap of the selected variable vs (phi, s).
   Rows with bnd == 0 (wall/limiter) are masked for the target-flux
   columns (angle, nv, nvT_gam, KparT) because those are not computed
   there. Use --rz for a poloidal R-Z view of one toroidal plane.

2. target_strike_profile.dat  -- per-plate profiles with signed distance
   from the strike point (sdist < 0: private flux region, sdist > 0:
   scrape-off layer).
   Columns: sdist  R  Z  angle  density  Ti  Te  Vpar  nv  nvT_gam
            KparT  psi  sputter
   Default: selected variable vs sdist, one line per plate.

Usage:
    python plot_target.py <file> [options]

Examples:
    # 2D surface map of parallel heat flux (phi vs s)
    python plot_target.py target_surface2d.dat -y KparT

    # 2D surface map of electron temperature, save to png
    python plot_target.py target_surface2d.dat -y Te -o te_surface.png

    # R-Z view of the phi = 90 deg plane
    python plot_target.py target_surface2d.dat -y Te --rz --phi 90

    # 1D density profile along targets vs strike distance
    python plot_target.py target_strike_profile.dat -y density

    # Ti vs sdist, outer target only (plate 1)
    python plot_target.py target_strike_profile.dat -y Ti --plate 1

    # List known variables
    python plot_target.py target_surface2d.dat --list-vars
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse
import sys

# ---------------------------------------------------------------------------
# Variable metadata: (description, unit_label, conversion_factor, display_name)
# ---------------------------------------------------------------------------
VARIABLE_INFO = {
    's':       ('Arc length along boundary',      'm',                      1.0,  r'$s$'),
    'strike_dist': ('Signed distance from strike pt', 'm',                 1.0,  r'$s_{\mathrm{strike}}$'),
    'R':       ('Major radius',                   'm',                      1.0,  r'$R$'),
    'Z':       ('Vertical position',              'm',                      1.0,  r'$Z$'),
    'psi':     ('Poloidal flux',                  'Wb/rad',                 1.0,  r'$\psi$'),
    'Psi_N':   ('Normalized poloidal flux',       '',                       1.0,  r'$\Psi_N$'),
    'rho':     ('Effective radius sqrt(Psi_N)',   '',                       1.0,  r'$\rho$'),
    'density': ('Mass density',                   r'$10^{20}\,\mathrm{m}^{-3}$', 1.0, r'$n_i$'),
    'Ti':      ('Ion temperature',                'keV',                    1.0,  r'$T_i$'),
    'Te':      ('Electron temperature',           'keV',                    1.0,  r'$T_e$'),
    'Vpar':    ('Parallel velocity',              'm/s',                    1.0,  r'$v_\parallel$'),
    'angle':   ('B-field incidence angle',        'deg',                    1.0,  r'$\theta$'),
    'nv':      ('Particle flux',                  r'$10^{20}\,\mathrm{m}^{-2}\mathrm{s}^{-1}$', 1.0, r'$nv_\parallel$'),
    'nvT_gam': ('Sheath heat flux',               r'$\mathrm{W}/\mathrm{m}^2$', 1.0, r'$\gamma n T v_\parallel$'),
    'KparT':   ('Parallel conduction flux',       r'$\mathrm{W}/\mathrm{m}^2$', 1.0, r'$\kappa_\parallel\nabla_\parallel T$'),
    'sputter': ('Sputtering yield',               '',                       1.0,  'Y'),
    'bnd':     ('Target flag (1 = target)',       '',                       1.0,  'bnd'),
}

# Columns only meaningful on target plates (bnd == 1)
TARGET_ONLY_COLS = ('angle', 'nv', 'nvT_gam', 'KparT')


def get_var_info(var):
    return VARIABLE_INFO.get(var, (var, '', 1.0, var))


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_surface2d(filepath):
    """Parse a target_surface2d.dat file (midplane2d-style layout).

    Returns
    -------
    names : list of str
    columns : list of dict
        Each dict has 'col' (1-based index), 'phi_deg' (float), and
        'data' (2D ndarray: npts x n_expr).
    """
    names = []
    columns = []
    current = None

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('# Column'):
                if current is not None and len(current['data']) > 0:
                    columns.append(current)
                tokens = stripped.split()
                phi_deg = None
                if '=' in stripped:
                    phi_deg = float(stripped.split('=')[1])
                current = {'col': int(tokens[2]), 'phi_deg': phi_deg, 'data': []}
                continue
            if stripped.startswith('#'):
                if not names:
                    names = stripped.lstrip('#').split()
                continue
            if current is not None:
                current['data'].append([float(x) for x in stripped.split()])

    if current is not None and len(current['data']) > 0:
        columns.append(current)

    if not names or not columns:
        print(f"ERROR: could not parse 2D surface file '{filepath}'.")
        sys.exit(1)

    npts = len(columns[0]['data'])
    for c in columns:
        c['data'] = np.array(c['data'])
        if len(c['data']) != npts:
            print(f"ERROR: column {c['col']} has {len(c['data'])} rows, "
                  f"expected {npts}.")
            sys.exit(1)
        if c['phi_deg'] is None:
            c['phi_deg'] = 360.0 * (c['col'] - 1) / len(columns)

    return names, columns


def parse_strike_profile(filepath):
    """Parse a target_strike_profile.dat file.

    Returns
    -------
    names : list of str
    plates : list of dict
        Each dict has 'plate_id' (int), 'strike_R'/'strike_Z' (float or None),
        and 'data' (2D ndarray: npts x n_expr).
    """
    names = []
    plates = []
    current = None

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('# Plate'):
                if current is not None and len(current['data']) > 0:
                    plates.append(current)
                sr = sz = None
                parts = stripped.split('=')
                if len(parts) == 3:
                    sr = float(parts[1].split()[0])
                    sz = float(parts[2].split()[0])
                current = {'plate_id': int(stripped.split()[2]),
                           'strike_R': sr, 'strike_Z': sz, 'data': []}
                continue
            if stripped.startswith('#'):
                if not names:
                    names = stripped.lstrip('#').split()
                continue
            if current is not None:
                current['data'].append([float(x) for x in stripped.split()])

    if current is not None and len(current['data']) > 0:
        plates.append(current)

    if not names or not plates:
        print(f"ERROR: could not parse strike profile file '{filepath}'.")
        sys.exit(1)

    ncol = len(plates[0]['data'][0])
    for p in plates:
        p['data'] = np.array(p['data'])
        if p['data'].shape[1] != ncol:
            print(f"ERROR: plate {p['plate_id']} has {p['data'].shape[1]} "
                  f"columns, expected {ncol}.")
            sys.exit(1)

    return names, plates


def parse_target_profile(filepath):
    """Parse a raw target_profile file (nearest-neighbor walk order)."""
    names = []
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                if not names:
                    names = stripped.lstrip('#').split()
                continue
            data.append([float(x) for x in stripped.split()])

    if not names or not data:
        print(f"ERROR: could not parse target profile file '{filepath}'.")
        sys.exit(1)

    data = np.array(data)
    if data.shape[1] > len(names):
        names = names + ['sputter']
    names = names[:data.shape[1]]
    return names, data


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def setup_plot_style():
    plt.rcParams.update({
        'font.size': 28,
        'font.weight': 'bold',
        'axes.linewidth': 4,
        'lines.linewidth': 3,
        'xtick.labelsize': 24,
        'ytick.labelsize': 24,
        'axes.labelsize': 32,
        'legend.fontsize': 22,
        'figure.dpi': 100,
        'text.usetex': False,
        'mathtext.default': 'regular',
    })


def plot_surface2d(names, columns, y_var='density', col_idx=None,
                   rz=False, title=None, output_file=None,
                   interactive=True, vmin=None, vmax=None):
    """Plot a 2D surface map: phi vs s colormap, or R-Z view."""
    if y_var not in names:
        print(f"ERROR: variable '{y_var}' not found. Available: {names}")
        sys.exit(1)

    setup_plot_style()

    y_idx = names.index(y_var)
    y_data = np.array([c['data'][:, y_idx] for c in columns])  # n_plane x npts
    phi_deg = np.array([c['phi_deg'] for c in columns])
    s_vals = columns[0]['data'][:, 0]

    if y_var in TARGET_ONLY_COLS and 'bnd' in names:
        bnd = columns[0]['data'][:, names.index('bnd')]
        y_data[:, bnd == 0] = np.nan  # wall rows: not computed

    info = get_var_info(y_var)
    y_label = info[3]
    y_unit = info[1]
    label = f'{y_label} ({y_unit})' if y_unit else y_label

    if rz:
        if col_idx is None:
            col_idx = len(columns) // 2
        c = columns[col_idx]
        R = c['data'][:, names.index('R')]
        Z = c['data'][:, names.index('Z')]
        vals = c['data'][:, y_idx]
        if y_var in TARGET_ONLY_COLS and 'bnd' in names:
            vals = vals.copy()
            vals[c['data'][:, names.index('bnd')] == 0] = np.nan
        sc = plt.scatter(R, Z, c=vals, cmap='viridis', s=60,
                         vmin=vmin, vmax=vmax)
        plt.xlabel(r'$R$ (m)')
        plt.ylabel(r'$Z$ (m)')
        plt.axis('equal')
        if y_var in TARGET_ONLY_COLS and 'bnd' in names:
            plt.scatter(R[c['data'][:, names.index('bnd')] == 0],
                        Z[c['data'][:, names.index('bnd')] == 0],
                        s=10, c='gray', label='wall')
            plt.legend(loc='best')
        plt.colorbar(sc).set_label(label)
        if title is None:
            title = f'phi = {c["phi_deg"]:.1f} deg'
    else:
        fig, ax = plt.subplots(figsize=(12.6, 10.0))
        fig.subplots_adjust(left=0.12, bottom=0.12, right=0.94, top=0.94)
        pc = ax.pcolormesh(phi_deg, s_vals, y_data.T,
                           shading='auto', cmap='viridis',
                           vmin=vmin, vmax=vmax)
        ax.set_xlabel(r'$\varphi$ (deg)')
        ax.set_ylabel(r'$s$ (m)')
        cb = plt.colorbar(pc, ax=ax)
        cb.set_label(label)
        if title:
            ax.set_title(title)

    if output_file:
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_file}")
    if interactive:
        plt.show()


def plot_strike_profile(names, plates, x_var='strike_dist', y_var='density',
                        plate_ids=None, title=None, output_file=None,
                        interactive=True):
    """Plot 1D profiles from target_strike_profile.dat, one line per plate."""
    if x_var not in names:
        print(f"ERROR: x-variable '{x_var}' not found. Available: {names}")
        sys.exit(1)
    if y_var not in names:
        print(f"ERROR: y-variable '{y_var}' not found. Available: {names}")
        sys.exit(1)

    setup_plot_style()

    x_idx = names.index(x_var)
    y_idx = names.index(y_var)

    fig, ax = plt.subplots(figsize=(12.6, 10.0))
    fig.subplots_adjust(left=0.12, bottom=0.12, right=0.94, top=0.94)

    plotted = 0
    for p in plates:
        if plate_ids is not None and p['plate_id'] not in plate_ids:
            continue
        x = p['data'][:, x_idx]
        y = p['data'][:, y_idx]
        # sort by x for a clean line
        order = np.argsort(x)
        ax.plot(x[order], y[order],
                label=f"plate {p['plate_id']}", linewidth=3)
        plotted += 1

    if plotted == 0:
        print("ERROR: no plates selected.")
        sys.exit(1)

    x_info = get_var_info(x_var)
    y_info = get_var_info(y_var)
    ax.set_xlabel(f'{x_info[3]} ({x_info[1]})' if x_info[1] else x_info[3])
    ax.set_ylabel(f'{y_info[3]} ({y_info[1]})' if y_info[1] else y_info[3])

    if title:
        ax.set_title(title)

    if plotted > 1:
        ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    if output_file:
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_file}")
    if interactive:
        plt.show()


def plot_target_profile(names, data, y_var='density', title=None,
                        output_file=None, interactive=True):
    """Plot raw target_profile (walk-ordered) vs arc index or R."""
    if y_var not in names:
        print(f"ERROR: y-variable '{y_var}' not found. Available: {names}")
        sys.exit(1)

    setup_plot_style()

    y_idx = names.index(y_var)
    x = np.arange(len(data))
    y = data[:, y_idx]

    info = get_var_info(y_var)
    fig, ax = plt.subplots(figsize=(12.6, 10.0))
    ax.plot(x, y, linewidth=3)
    ax.set_xlabel('point index along boundary')
    ax.set_ylabel(f'{info[3]} ({info[1]})' if info[1] else info[3])
    if title:
        ax.set_title(title)
    ax.grid(True, alpha=0.3)

    if output_file:
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_file}")
    if interactive:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Plot jorek2_target2vtkn ASCII outputs '
                    '(target_surface2d.dat, target_strike_profile.dat, '
                    'target_profile).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s target_surface2d.dat -y KparT
  %(prog)s target_surface2d.dat -y Te --rz --phi 90
  %(prog)s target_strike_profile.dat -y density
  %(prog)s target_strike_profile.dat -y Ti --plate 1
        """)
    parser.add_argument('file', nargs='?',
                        help='Path to a target2vtkn output file')
    parser.add_argument('-x', '--xvar', default='strike_dist',
                        help='x-variable for 1D profiles (default: strike_dist)')
    parser.add_argument('-y', '--yvar', default=None,
                        help='variable to plot (default: density)')
    parser.add_argument('--plate', type=int, nargs='+',
                        help='1D mode: restrict to plate ID(s)')
    parser.add_argument('--rz', action='store_true',
                        help='2D mode: R-Z poloidal view instead of '
                             'phi-s colormap')
    parser.add_argument('--col', type=int,
                        help='2D mode: select toroidal column index '
                             '(1-based) for --rz')
    parser.add_argument('--phi', type=float, dest='phi_deg',
                        help='2D mode: select column closest to this '
                             'toroidal angle in degrees')
    parser.add_argument('--vmin', type=float, help='Color scale minimum')
    parser.add_argument('--vmax', type=float, help='Color scale maximum')
    parser.add_argument('-o', '--output', help='Save plot to file')
    parser.add_argument('--title', help='Plot title')
    parser.add_argument('--no-show', action='store_true',
                        help='Do not display interactive window')
    parser.add_argument('--agg', action='store_true',
                        help='Use non-interactive Agg backend')
    parser.add_argument('--list-vars', action='store_true',
                        help='List known variables and exit')
    args = parser.parse_args()

    if args.agg:
        matplotlib.use('Agg')

    if args.list_vars:
        print("\nKnown variables with display names and default units:")
        for v, info in VARIABLE_INFO.items():
            unit = info[1] if info[1] else 'dimensionless'
            print(f"  {v:12s} -> {info[3]:35s} [{unit}]")
        return

    if args.file is None:
        parser.print_usage()
        print("ERROR: input file is required (unless --list-vars is used).")
        sys.exit(1)

    with open(args.file, 'r') as f:
        head = ''.join(f.readline() for _ in range(8))
    is_2d = '# Column' in head

    if is_2d:
        names, columns = parse_surface2d(args.file)
        print(f"File: {args.file}  (2D surface, {len(columns)} toroidal "
              f"planes x {len(columns[0]['data'])} boundary points)")
        print(f"Variables: {names}")

        y_var = args.yvar if args.yvar else 'density'
        if y_var not in names:
            print(f"ERROR: variable '{y_var}' not found. Available: {names}")
            sys.exit(1)

        col_idx = None
        if args.col is not None:
            if not (1 <= args.col <= len(columns)):
                print(f"ERROR: --col {args.col} out of range "
                      f"(1..{len(columns)}).")
                sys.exit(1)
            col_idx = args.col - 1
        if args.phi_deg is not None:
            phis = [c['phi_deg'] for c in columns]
            col_idx = min(range(len(phis)),
                          key=lambda i: abs(phis[i] - args.phi_deg))
            print(f"phi = {args.phi_deg} deg -> column {col_idx + 1} "
                  f"(phi = {phis[col_idx]:.2f} deg)")

        plot_surface2d(names, columns, y_var=y_var, col_idx=col_idx,
                       rz=args.rz, title=args.title,
                       output_file=args.output,
                       interactive=not args.no_show,
                       vmin=args.vmin, vmax=args.vmax)
    else:
        names, plates = parse_strike_profile(args.file)
        print(f"File: {args.file}  (1D strike profiles, "
              f"{len(plates)} plates)")
        print(f"Variables: {names}")

        y_var = args.yvar if args.yvar else 'density'
        if y_var not in names:
            print(f"ERROR: variable '{y_var}' not found. Available: {names}")
            sys.exit(1)

        plot_strike_profile(names, plates, x_var=args.xvar, y_var=y_var,
                            plate_ids=args.plate, title=args.title,
                            output_file=args.output,
                            interactive=not args.no_show)


if __name__ == '__main__':
    main()
