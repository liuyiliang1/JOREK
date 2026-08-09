#!/usr/bin/env python3
"""
Plot JOREK2 post-processed midplane data.

Read exprs_midplane_*.dat files produced by jorek2_postproc and plot
any variable against any other, with appropriate unit conversions.

Default: T_e vs Psi_N on the outer midplane (low-field side, LFS).
The midplane profile goes inner wall -> magnetic axis -> outer wall.
By default only the low-field side (outer) is shown.

Usage:
    python plot_midplane.py <midplane_file> [options]

Examples:
    # T_e vs Psi_N, last step, LFS (default)
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat

    # ne vs R, all steps, both sides
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat -x R -y ne --all --side both

    # currdens vs Psi_N, single step
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat -y currdens --step 500

    # T_e vs rho = sqrt(Psi_N) (effective radius, linear in minor radius)
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat -x rho

    # List available variables
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat --list-vars
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse
import sys
import os

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
MU_ZERO = 4.0 * np.pi * 1.0e-7
EL_CHG  = 1.602176634e-19
AMU     = 1.66053906660e-27

# ---------------------------------------------------------------------------
# Variable metadata: (description, unit_label, conversion_factor, display_name)
# ---------------------------------------------------------------------------
VARIABLE_INFO = {
    'R':              ('Major radius',                   'm',                     1.0,  r'$R$'),
    'ne':             ('Electron density',                r'$10^{19}\,\mathrm{m}^{-3}$', 1e-19, r'$n_e$'),
    'Psi_N':          ('Normalized poloidal flux',        '',                      1.0,  r'$\Psi_N$'),
    'T_i':            ('Ion temperature',                 'keV',                   1e-3, r'$T_i$'),
    't':              ('Time',                            's',                     1.0,  r'$t$'),
    'pres':           ('Total pressure',                  'kPa',                   1e-3, r'$P$'),
    'currdens':       ('Toroidal current density',        r'$\mathrm{MA}/\mathrm{m}^2$', 1e-6, r'$J_\phi$'),
    'Er':             ('Radial electric field',           r'$\mathrm{kV}/\mathrm{m}$', 1e-3, r'$E_r$'),
    'gradPdotCurv':   (r'$\nabla P \cdot \kappa$',        r'$\mathrm{N}/\mathrm{m}^3$', 1.0,  r'$\nabla P\cdot\kappa$'),
    'dprof':          ('Particle diffusivity',            r'$\mathrm{m}^2/\mathrm{s}$', 1.0,  r'$D_\perp$'),
    'zkiprof':        ('Ion heat diffusivity',            r'$\mathrm{m}^2/\mathrm{s}$', 1.0,  r'$\kappa_\perp$'),
    'ki_neo':         ('Neoclassical heat diffusivity',   r'$\mathrm{m}^2/\mathrm{s}$', 1.0,  r'$\kappa_{\mathrm{neo}}$'),
    'mu_neo':         ('Neoclassical friction coeff.',    '',                      1.0,  r'$\mu_{\mathrm{neo}}$'),
    'T_e':            ('Electron temperature',            'keV',                   1e-3, r'$T_e$'),
    'eta_T':          ('Resistivity',                     r'$\Omega\cdot\mathrm{m}$', 1.0,  r'$\eta$'),
    'visco_T':        ('Viscosity',                       r'$\mathrm{m}^2/\mathrm{s}$', 1.0,  r'$\nu$'),
    'J_bootstrap':    ('Bootstrap current density',       r'$\mathrm{MA}/\mathrm{m}^2$', 1e-6, r'$J_{\mathrm{bs}}$'),
    'JxB_R':          (r'$\mathbf{J}\times\mathbf{B}|_R$', r'$\mathrm{N}/\mathrm{m}^3$', 1.0,  r'$\mathbf{J}\times\mathbf{B}|_R$'),
    'gradP_R':        (r'$\nabla P|_R$',                  r'$\mathrm{N}/\mathrm{m}^3$', 1.0,  r'$\nabla P|_R$'),
    'vpar':           ('Parallel velocity',               r'$\mathrm{m}/\mathrm{s}$', 1.0,  r'$v_\parallel$'),
    'V_ExB_pol':      ('Poloidal ExB velocity',           r'$\mathrm{m}/\mathrm{s}$', 1.0,  r'$V_{\mathrm{ExB},\mathrm{pol}}$'),
    'V_ExB_R':        ('ExB velocity R',                  r'$\mathrm{m}/\mathrm{s}$', 1.0,  r'$V_{\mathrm{ExB},R}$'),
    'V_ExB_Z':        ('ExB velocity Z',                  r'$\mathrm{m}/\mathrm{s}$', 1.0,  r'$V_{\mathrm{ExB},Z}$'),
}

# ---------------------------------------------------------------------------
# File reader
# ---------------------------------------------------------------------------

def parse_midplane_file(filepath):
    """Parse a jorek2_postproc midplane output file.

    Returns
    -------
    var_names : list of str
        Column names from the header.
    time_steps : list of dict
        Each dict has keys 'step' (int), 'data' (2D ndarray: npts x nvars).
    """
    var_names = []
    time_steps = []
    current_step = None
    current_data = []

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                if current_step is not None and current_data:
                    time_steps.append({
                        'step': current_step,
                        'data': np.array(current_data)
                    })
                    current_data = []
                    current_step = None
                continue

            if stripped.startswith('# time step'):
                if current_step is not None and current_data:
                    time_steps.append({
                        'step': current_step,
                        'data': np.array(current_data)
                    })
                    current_data = []
                parts = stripped.split('#')
                current_step = int(parts[-1].strip())
                continue

            if stripped.startswith('#'):
                if not var_names:
                    var_names = stripped.lstrip('#').split()
                continue

            if current_step is not None:
                values = [float(x) for x in stripped.split()]
                current_data.append(values)

    if current_step is not None and current_data:
        time_steps.append({
            'step': current_step,
            'data': np.array(current_data)
        })

    return var_names, time_steps


# ---------------------------------------------------------------------------
# Midplane side splitting
# ---------------------------------------------------------------------------

def split_outer_midplane(data, psi_idx):
    """Extract low-field side (outer midplane): magnetic axis -> outer wall."""
    psi = data[:, psi_idx]
    idx_min = np.argmin(psi)
    return data[idx_min:, :]


def split_inner_midplane(data, psi_idx):
    """Extract high-field side (inner midplane): inner wall -> magnetic axis."""
    psi = data[:, psi_idx]
    idx_min = np.argmin(psi)
    return data[:idx_min + 1, :]


# ---------------------------------------------------------------------------
# Unit conversion
# ---------------------------------------------------------------------------

def compute_jorek_conversion_factors(central_density=0.973, central_mass=2.0):
    """Compute conversion factors from JOREK normalized to SI units."""
    n0_SI = central_density * 1.0e20
    rho_norm = n0_SI * central_mass * AMU
    fact_time = np.sqrt(MU_ZERO * rho_norm)
    fact_T = 1.0 / (MU_ZERO * n0_SI * EL_CHG)
    fact_resistiv = np.sqrt(MU_ZERO / rho_norm)
    fact_rho = n0_SI * central_mass * AMU

    return {
        'fact_time':     fact_time,
        'fact_T':        fact_T,
        'fact_ne':       n0_SI,
        'fact_pres':     1.0 / MU_ZERO,
        'fact_currdens': 1.0 / MU_ZERO,
        'fact_dprof':    1.0 / fact_time,
        'fact_zkiprof':  1.0 / (fact_resistiv * (5.0/3.0 - 1.0) * fact_rho),
        'fact_eta':      fact_resistiv,
        'fact_visco':    1.0 / fact_resistiv,
        'fact_vpar':     1.0 / fact_time,
    }


def convert_variable(data, var_name, conversion_factor=None,
                     central_density=0.973, central_mass=2.0,
                     assume_si=True):
    """Apply unit conversion to a variable's data."""
    info = VARIABLE_INFO.get(var_name, (var_name, '', 1.0, var_name))

    if conversion_factor is not None:
        return data * conversion_factor, info[1]

    if assume_si:
        return data * info[2], info[1]
    else:
        factors = compute_jorek_conversion_factors(central_density, central_mass)

        if var_name in ('T_e', 'T_i'):
            si = data * factors['fact_T']
            display = si * 1e-3
            unit_label = 'keV'
        elif var_name == 'ne':
            si = data * factors['fact_ne']
            display = si * 1e-19
            unit_label = r'$10^{19}\,\mathrm{m}^{-3}$'
        elif var_name == 'pres':
            si = data * factors['fact_pres']
            display = si * 1e-3
            unit_label = 'kPa'
        elif var_name == 'currdens':
            si = data * factors['fact_currdens']
            display = si * 1e-6
            unit_label = r'$\mathrm{MA}/\mathrm{m}^2$'
        elif var_name in ('dprof',):
            si = data * factors['fact_dprof']
            display = si
            unit_label = r'$\mathrm{m}^2/\mathrm{s}$'
        elif var_name in ('zkiprof', 'ki_neo'):
            si = data * factors['fact_zkiprof']
            display = si
            unit_label = r'$\mathrm{m}^2/\mathrm{s}$'
        elif var_name == 't':
            display = data * factors['fact_time']
            unit_label = 's'
        elif var_name == 'vpar':
            display = data * factors['fact_vpar']
            unit_label = 'm/s'
        else:
            display = data * info[2]
            unit_label = info[1]

        return display, unit_label


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def setup_plot_style():
    """Configure matplotlib to match MATLAB reference style."""
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


def plot_midplane(var_names, time_steps, x_var='Psi_N', y_var='T_e',
                  selected_steps=None, conversion_factor_x=None,
                  conversion_factor_y=None, assume_si=True,
                  central_density=0.973, central_mass=2.0,
                  side='outer',
                  title=None, output_file=None, interactive=True):
    """Plot selected variables from midplane data."""
    setup_plot_style()

    # 'rho' is a derived x-variable: rho = sqrt(Psi_N)
    is_rho = (x_var == 'rho')

    if is_rho:
        if 'Psi_N' not in var_names:
            print("ERROR: x-variable 'rho' requires a 'Psi_N' column in the file.")
            sys.exit(1)
        if 'rho' in var_names:
            # JOREK postproc names the mass-density column 'rho' as well;
            # as an x-variable we always mean the derived radial coordinate.
            print("NOTE: file also contains a 'rho' column (JOREK mass density); "
                  "'-x rho' uses derived sqrt(Psi_N). Density is still "
                  "available as a y-variable.")
        x_idx = var_names.index('Psi_N')
    else:
        if x_var not in var_names:
            print(f"ERROR: x-variable '{x_var}' not found. Available: {var_names}")
            sys.exit(1)
        x_idx = var_names.index(x_var)

    if y_var not in var_names:
        print(f"ERROR: y-variable '{y_var}' not found. Available: {var_names}")
        sys.exit(1)

    y_idx = var_names.index(y_var)

    # Determine Psi_N index for side filtering
    psi_idx = None
    if side in ('outer', 'inner') and 'Psi_N' in var_names:
        psi_idx = var_names.index('Psi_N')

    fig, ax = plt.subplots(figsize=(12.6, 10.0))
    fig.subplots_adjust(left=0.12, bottom=0.12, right=0.94, top=0.94)
    ax.set_box_aspect(1.15 / 1.45)

    steps_to_plot = time_steps
    if selected_steps is not None:
        steps_to_plot = [ts for ts in time_steps if ts['step'] in selected_steps]

    if not steps_to_plot:
        print("ERROR: No matching time steps found.")
        sys.exit(1)

    cmap = plt.cm.viridis
    n_steps = len(steps_to_plot)
    for i, ts in enumerate(steps_to_plot):
        data = ts['data']

        if side == 'outer' and psi_idx is not None:
            data = split_outer_midplane(data, psi_idx)
        elif side == 'inner' and psi_idx is not None:
            data = split_inner_midplane(data, psi_idx)

        x_raw = data[:, x_idx]
        y_raw = data[:, y_idx]
        if is_rho:
            # rho = sqrt(Psi_N): clip at 0 to guard against tiny negative
            # values from numerical noise near the magnetic axis.
            x_data = np.sqrt(np.clip(x_raw, 0.0, None))
            x_unit = ''
        else:
            x_data, x_unit = convert_variable(
                x_raw, x_var, conversion_factor_x, central_density, central_mass, assume_si)
        y_data, y_unit = convert_variable(
            y_raw, y_var, conversion_factor_y, central_density, central_mass, assume_si)

        color = cmap(float(i) / max(n_steps - 1, 1))
        side_label = {'outer': ' (LFS)', 'inner': ' (HFS)', 'both': ''}
        label = f"step {ts['step']}{side_label.get(side, '')}"
        ax.plot(x_data, y_data, color=color, label=label, linewidth=3)

    if is_rho:
        x_label = r'$\rho$'
    else:
        x_label = VARIABLE_INFO.get(x_var, (x_var, '', 1.0, x_var))[3]
    y_label = VARIABLE_INFO.get(y_var, (y_var, '', 1.0, y_var))[3]

    if x_unit:
        ax.set_xlabel(f'{x_label} ({x_unit})')
    else:
        ax.set_xlabel(x_label)

    if y_unit:
        ax.set_ylabel(f'{y_label} ({y_unit})')
    else:
        ax.set_ylabel(y_label)

    if title:
        ax.set_title(title)

    if n_steps > 1 and n_steps <= 20:
        ax.legend(loc='best')
    elif n_steps > 20:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(
            vmin=steps_to_plot[0]['step'], vmax=steps_to_plot[-1]['step']))
        plt.colorbar(sm, ax=ax).set_label('Time step')

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
        description='Plot JOREK2 midplane data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s exprs_midplane_s000000..000700.dat
  %(prog)s exprs_midplane_s000000..000700.dat -x R -y ne --all --side both
  %(prog)s exprs_midplane_s000000..000700.dat -y currdens --step 700
  %(prog)s exprs_midplane_s000000..000700.dat -x rho
  %(prog)s exprs_midplane_s000000..000700.dat -o plot.png
  %(prog)s exprs_midplane_s000000..000700.dat --jorek-units --n0 1.0 --mass 2.0
        """)
    parser.add_argument('file', help='Path to midplane .dat file')
    parser.add_argument('-x', '--xvar', default='Psi_N',
                        help='Variable for x-axis (default: Psi_N); "rho" plots '
                             'sqrt(Psi_N), the normalized effective radius')
    parser.add_argument('-y', '--yvar', default='T_e',
                        help='Variable for y-axis (default: T_e)')
    parser.add_argument('--step', type=int, nargs='+',
                        help='Plot specific time step(s)')
    parser.add_argument('--step-index', type=int, nargs='+',
                        help='Plot specific step index/indices (0-based)')
    parser.add_argument('--all', action='store_true',
                        help='Plot all time steps (default: last step only)')
    parser.add_argument('--side', choices=['outer', 'inner', 'both'], default='outer',
                        help='Midplane side: outer/LFS (default), inner/HFS, or both')
    parser.add_argument('-o', '--output', help='Save plot to file')
    parser.add_argument('--title', help='Plot title')
    parser.add_argument('--no-show', action='store_true',
                        help='Do not display interactive window')
    parser.add_argument('--agg', action='store_true',
                        help='Use non-interactive Agg backend')

    unit_group = parser.add_argument_group('Unit conversion')
    unit_group.add_argument('--assume-si', action='store_true', default=True,
                            help='Assume raw data is in SI units [default]')
    unit_group.add_argument('--jorek-units', action='store_true',
                            help='Raw data is in JOREK normalized units')
    unit_group.add_argument('--n0', type=float, default=0.973,
                            help='Central density (in 10^20 m^-3)')
    unit_group.add_argument('--mass', type=float, default=2.0,
                            help='Average ion mass in AMU')
    unit_group.add_argument('--cx', type=float,
                            help='Direct conversion factor for x-variable')
    unit_group.add_argument('--cy', type=float,
                            help='Direct conversion factor for y-variable')
    parser.add_argument('--list-vars', action='store_true',
                        help='List available variables in the file and exit')

    args = parser.parse_args()

    if args.agg:
        matplotlib.use('Agg')

    var_names, time_steps = parse_midplane_file(args.file)

    if not var_names:
        print("ERROR: Could not parse variable names from header.")
        sys.exit(1)
    if not time_steps:
        print("ERROR: No time step data found in file.")
        sys.exit(1)

    print(f"File: {args.file}")
    print(f"Variables: {var_names}")
    print(f"Time steps: {len(time_steps)} steps "
          f"({time_steps[0]['step']} to {time_steps[-1]['step']})")
    print(f"Points per step: {time_steps[0]['data'].shape[0]}")

    if args.list_vars:
        print("\nAvailable variables with display names and default units:")
        for v in var_names:
            info = VARIABLE_INFO.get(v, (v, '?', 1.0, v))
            unit = info[1] if info[1] else 'dimensionless'
            print(f"  {v:20s} -> {info[3]:30s} [{unit}]")
        print(f"  {'rho':20s} -> {'sqrt(Psi_N), effective radius':30s} [derived]")
        return

    assume_si = not args.jorek_units

    selected_steps = None
    if args.step is not None:
        selected_steps = args.step
    elif args.step_index is not None:
        selected_steps = [time_steps[i]['step'] for i in args.step_index
                          if 0 <= i < len(time_steps)]
    elif not args.all:
        selected_steps = [time_steps[-1]['step']]
        print(f"Plotting last step only (use --all for all {len(time_steps)} steps)")

    plot_midplane(
        var_names, time_steps,
        x_var=args.xvar, y_var=args.yvar,
        selected_steps=selected_steps,
        conversion_factor_x=args.cx,
        conversion_factor_y=args.cy,
        assume_si=assume_si,
        central_density=args.n0,
        central_mass=args.mass,
        side=args.side,
        title=args.title,
        output_file=args.output,
        interactive=not args.no_show,
    )


if __name__ == '__main__':
    main()
