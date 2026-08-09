#!/usr/bin/env python3
"""
Plot JOREK2 post-processed midplane data.

Read the exprs_midplane_*.dat files produced by jorek2_postproc and plot
any variable against any other, with appropriate unit conversions.

Usage:
    python plot_midplane.py <midplane_file> [options]

Examples:
    # Plot T_e vs Psi_N for the last time step (default)
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat

    # Plot ne vs R for all time steps
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat -x R -y ne --all

    # Plot currdens vs Psi_N for a single time step
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat -y currdens --step 500

    # Plot specific steps
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat --step 0 350 700

    # Save to file instead of showing
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat -o my_plot.png

    # Custom unit reference parameters (for JOREK-unit data)
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat --jorek-units --n0 1.0 --mass 2.0
"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Can be overridden with --agg for headless
import matplotlib.pyplot as plt
import argparse
import sys
import os

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
MU_ZERO = 4.0 * np.pi * 1.0e-7          # Vacuum permeability [H/m]
EL_CHG  = 1.602176634e-19                # Elementary charge [C]
AMU     = 1.66053906660e-27              # Atomic mass unit [kg]

# ---------------------------------------------------------------------------
# Variable metadata: (description, unit_label, conversion_factor, display_name)
#   display_name: LaTeX-formatted name for axis labels
#   conversion_factor: multiply raw data by this to get display units
# ---------------------------------------------------------------------------
VARIABLE_INFO = {
    'R':              ('Major radius',                'm',                    1.0,  '$R$'),
    'ne':             ('Electron density',             '$10^{19}\\,m^{-3}$',  1e-19,'$n_e$'),
    'Psi_N':          ('Normalized poloidal flux',     '',                     1.0,  '$\Psi_N$'),
    'T_i':            ('Ion temperature',              'keV',                  1e-3, '$T_i$'),
    't':              ('Time',                         's',                    1.0,  '$t$'),
    'pres':           ('Total pressure',               'kPa',                  1e-3, '$P$'),
    'currdens':       ('Toroidal current density',     '$MA/m^2$',             1e-6, '$J_\phi$'),
    'Er':             ('Radial electric field',        '$kV/m$',               1e-3, '$E_r$'),
    'gradPdotCurv':   ('grad P dot curvature',         '$N/m^3$',              1.0,  r'$\nabla P\cdot\kappa$'),
    'dprof':          ('Particle diffusivity',         '$m^2/s$',              1.0,  r'$D_\perp$'),
    'zkiprof':        ('Ion heat diffusivity',         '$m^2/s$',              1.0,  r'$\kappa_\perp$'),
    'ki_neo':         ('Neoclassical heat diffusivity','$m^2/s$',              1.0,  r'$\kappa_{neo}$'),
    'mu_neo':         ('Neoclassical friction coeff.', '',                     1.0,  r'$\mu_{neo}$'),
    'T_e':            ('Electron temperature',         'keV',                  1e-3, r'$T_e$'),
    'eta_T':          ('Resistivity',                  r'$\Omega\cdot m$',     1.0,  r'$\eta$'),
    'visco_T':        ('Viscosity',                    '$m^2/s$',              1.0,  r'$\nu$'),
    'J_bootstrap':    ('Bootstrap current density',    '$MA/m^2$',             1e-6, r'$J_{bs}$'),
    'JxB_R':          ('JxB force (R component)',      '$N/m^3$',              1.0,  r'$J\times B|_R$'),
    'gradP_R':        ('Pressure gradient (R comp.)',  '$N/m^3$',              1.0,  r'$\nabla P|_R$'),
    'vpar':           ('Parallel velocity',            '$m/s$',                1.0,  r'$v_\parallel$'),
    'V_ExB_pol':      ('Poloidal ExB velocity',        '$m/s$',                1.0,  '$V_{ExB,pol}$'),
    'V_ExB_R':        ('ExB velocity (R component)',   '$m/s$',                1.0,  '$V_{ExB,R}$'),
    'V_ExB_Z':        ('ExB velocity (Z component)',   '$m/s$',                1.0,  '$V_{ExB,Z}$'),
}

# Variables that may be in JOREK units and need physics-based conversion
# If the file was generated with "set units 1" (JOREK units), the raw data is
# in normalized units and needs conversion factors derived from reference
# parameters: n0, m_i, etc.
JOREK_UNIT_VARIABLES = {
    # variable: (scale_with_n0, scale_with_T0_factor)
    # T0_factor = 1/(MU_ZERO * n0 * 1e20 * EL_CHG) converts JOREK T to eV
    'T_e':  ('temperature',),
    'T_i':  ('temperature',),
    'ne':   ('density',),
    'pres': ('pressure',),
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
                # Blank line: flush current time step if data exists
                if current_step is not None and current_data:
                    time_steps.append({
                        'step': current_step,
                        'data': np.array(current_data)
                    })
                    current_data = []
                    current_step = None
                continue

            if stripped.startswith('# time step'):
                # Flush previous step if any
                if current_step is not None and current_data:
                    time_steps.append({
                        'step': current_step,
                        'data': np.array(current_data)
                    })
                    current_data = []
                # Parse step number: "# time step #000100"
                parts = stripped.split('#')
                current_step = int(parts[-1].strip())
                continue

            if stripped.startswith('#'):
                # Header line: variable names
                if not var_names:
                    var_names = stripped.lstrip('#').split()
                continue

            # Data line
            if current_step is not None:
                values = [float(x) for x in stripped.split()]
                current_data.append(values)

    # Don't forget the last step
    if current_step is not None and current_data:
        time_steps.append({
            'step': current_step,
            'data': np.array(current_data)
        })

    return var_names, time_steps


# ---------------------------------------------------------------------------
# Unit conversion helpers
# ---------------------------------------------------------------------------

def compute_jorek_conversion_factors(central_density=0.973, central_mass=2.0):
    """Compute conversion factors from JOREK normalized to SI units.

    Parameters
    ----------
    central_density : float
        n0 in units of 10^20 m^-3 (from JOREK namelist).
    central_mass : float
        Average ion mass in AMU (from JOREK namelist).

    Returns
    -------
    dict of conversion factors
    """
    n0_SI = central_density * 1.0e20                         # [m^-3]
    rho_norm = n0_SI * central_mass * AMU                    # [kg/m^3]
    fact_time = np.sqrt(MU_ZERO * rho_norm)                  # [s]
    fact_T = 1.0 / (MU_ZERO * n0_SI * EL_CHG)                # T [JOREK] -> T [eV]
    fact_resistiv = np.sqrt(MU_ZERO / rho_norm)
    fact_rho = n0_SI * central_mass * AMU

    return {
        'fact_time':     fact_time,       # time: JOREK -> s
        'fact_T':        fact_T,          # Te/Ti: JOREK -> eV
        'fact_ne':       n0_SI,           # ne: JOREK -> m^-3
        'fact_pres':     1.0 / MU_ZERO,   # pres: JOREK -> Pa
        'fact_currdens': 1.0 / MU_ZERO,   # currdens: JOREK -> ... (approximate)
        'fact_dprof':    1.0 / fact_time, # dprof: JOREK -> m^2/s
        'fact_zkiprof':  1.0 / (fact_resistiv * (5.0/3.0 - 1.0) * fact_rho),
        'fact_eta':      fact_resistiv,
        'fact_visco':    1.0 / fact_resistiv,
        'fact_vpar':     1.0 / fact_time, # approximate
    }


def convert_variable(data, var_name, conversion_factor=None,
                     central_density=0.973, central_mass=2.0,
                     assume_si=True):
    """Apply unit conversion to a variable's data.

    Parameters
    ----------
    data : ndarray
        Raw data values.
    var_name : str
        Variable name.
    conversion_factor : float or None
        If given, directly multiply data by this factor.
        If None, use the default from VARIABLE_INFO.
    assume_si : bool
        If True, assume the data is already in SI and only scale to display
        units (e.g. eV -> keV, m^-3 -> 10^19 m^-3).
        If False, convert from JOREK units to display units using physics
        reference parameters.

    Returns
    -------
    converted_data : ndarray
    unit_label : str
    """
    info = VARIABLE_INFO.get(var_name, (var_name, '', 1.0, var_name))

    if conversion_factor is not None:
        return data * conversion_factor, info[1]

    if assume_si:
        # Data is already in SI units, just scale to display units
        return data * info[2], info[1]
    else:
        # Data is in JOREK normalized units, convert to SI first
        factors = compute_jorek_conversion_factors(central_density, central_mass)

        # Map variable names to conversion chains
        if var_name in ('T_e', 'T_i'):
            # JOREK -> eV -> keV
            si = data * factors['fact_T']        # eV
            display = si * 1e-3                  # keV
            unit_label = 'keV'
        elif var_name == 'ne':
            si = data * factors['fact_ne']        # m^-3
            display = si * 1e-19                  # 10^19 m^-3
            unit_label = '10^{19} m^{-3}'
        elif var_name == 'pres':
            si = data * factors['fact_pres']      # Pa
            display = si * 1e-3                   # kPa
            unit_label = 'kPa'
        elif var_name == 'currdens':
            si = data * factors['fact_currdens']  # A/m^2
            display = si * 1e-6                   # MA/m^2
            unit_label = 'MA/m^2'
        elif var_name in ('dprof',):
            si = data * factors['fact_dprof']     # m^2/s
            display = si
            unit_label = 'm^2/s'
        elif var_name in ('zkiprof', 'ki_neo'):
            si = data * factors['fact_zkiprof']   # m^2/s
            display = si
            unit_label = 'm^2/s'
        elif var_name == 'eta_T':
            si = data * factors['fact_eta']
            display = si
            unit_label = 'Ohm*m'
        elif var_name == 'visco_T':
            si = data * factors['fact_visco']
            display = si
            unit_label = 'm^2/s'
        elif var_name == 't':
            display = data * factors['fact_time']
            unit_label = 's'
        elif var_name == 'vpar':
            display = data * factors['fact_vpar']
            unit_label = 'm/s'
        elif var_name == 'Er':
            si = data * factors['fact_T']  # rough; Er uses F0/fact_time
            display = si
            unit_label = 'V/m'
        else:
            display = data * info[2]
            unit_label = info[1]

        return display, unit_label


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def setup_plot_style():
    """Configure matplotlib to match the reference MATLAB style."""
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
    })


def split_outer_midplane(data, psi_idx):
    """Extract the low-field side (outer midplane) from full midplane data.

    The midplane profile goes: inner wall -> magnetic axis (Psi_N min) -> outer wall.
    Returns data from the magnetic axis to the outer wall.
    """
    psi = data[:, psi_idx]
    idx_min = np.argmin(psi)
    return data[idx_min:, :]


def split_inner_midplane(data, psi_idx):
    """Extract the high-field side (inner midplane) from full midplane data.

    Returns data from the inner wall to the magnetic axis.
    """
    psi = data[:, psi_idx]
    idx_min = np.argmin(psi)
    return data[:idx_min + 1, :]


def plot_midplane(var_names, time_steps, x_var='R', y_var='T_e',
                  selected_steps=None, conversion_factor_x=None,
                  conversion_factor_y=None, assume_si=True,
                  central_density=0.973, central_mass=2.0,
                  side='outer',
                  title=None, output_file=None, interactive=True):
    """Plot selected variables from midplane data.

    Parameters
    ----------
    var_names : list of str
        Column names.
    time_steps : list of dict
        Parsed time step data.
    x_var, y_var : str
        Variable names for x and y axes.
    selected_steps : list of int or None
        Which time steps to plot. None = all steps.
    conversion_factor_x, conversion_factor_y : float or None
        Direct scaling factors.
    assume_si : bool
        If True, assume raw data is in SI.
    central_density, central_mass : float
        Reference parameters for JOREK -> SI conversion.
    side : str
        Which side of the midplane: 'outer' (low-field, default), 'inner' (high-field),
        or 'both' (full midplane).
    title : str or None
        Plot title.
    output_file : str or None
        If set, save to this file instead of showing.
    interactive : bool
        If True, show the plot window.
    """
    setup_plot_style()

    if x_var not in var_names:
        print(f"ERROR: x-variable '{x_var}' not found. Available: {var_names}")
        sys.exit(1)
    if y_var not in var_names:
        print(f"ERROR: y-variable '{y_var}' not found. Available: {var_names}")
        sys.exit(1)

    x_idx = var_names.index(x_var)
    y_idx = var_names.index(y_var)

    # Create figure with size matching MATLAB default
    fig, ax = plt.subplots(figsize=(12.6, 10.0))

    # Configure axes geometry (matching MATLAB settings)
    # ax.OuterPosition = [0.06, 0.06, 1, 1]
    fig.subplots_adjust(left=0.12, bottom=0.12, right=0.94, top=0.94)
    # ax.PlotBoxAspectRatio = [1.45, 1.15, 1]
    ax.set_box_aspect(1.15/1.45)

    # Plot each time step
    steps_to_plot = time_steps
    if selected_steps is not None:
        steps_to_plot = [ts for ts in time_steps if ts['step'] in selected_steps]

    if not steps_to_plot:
        print("ERROR: No matching time steps found.")
        sys.exit(1)

    # Determine Psi_N index for side filtering
    psi_idx = None
    if side in ('outer', 'inner') and 'Psi_N' in var_names:
        psi_idx = var_names.index('Psi_N')

    cmap = plt.cm.viridis
    n_steps = len(steps_to_plot)
    for i, ts in enumerate(steps_to_plot):
        data = ts['data']

        # Filter by midplane side (low-field / high-field)
        if side == 'outer' and psi_idx is not None:
            data = split_outer_midplane(data, psi_idx)
        elif side == 'inner' and psi_idx is not None:
            data = split_inner_midplane(data, psi_idx)

        x_raw = data[:, x_idx]
        y_raw = data[:, y_idx]

        # Apply conversions
        x_data, x_unit = convert_variable(
            x_raw, x_var, conversion_factor_x, central_density, central_mass, assume_si)
        y_data, y_unit = convert_variable(
            y_raw, y_var, conversion_factor_y, central_density, central_mass, assume_si)

        color = cmap(float(i) / max(n_steps - 1, 1))
        side_label = {None: '', 'outer': ' (LFS)', 'inner': ' (HFS)', 'both': ''}
        label = f"step {ts['step']}{side_label.get(side, '')}"
        ax.plot(x_data, y_data, color=color, label=label, linewidth=3)

    # Axis labels with LaTeX-formatted variable names and units
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

    # Legend if multiple steps
    if n_steps > 1 and n_steps <= 20:
        ax.legend(loc='best')
    elif n_steps > 20:
        # Add colorbar instead of legend for many steps
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(
            vmin=steps_to_plot[0]['step'], vmax=steps_to_plot[-1]['step']))
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label('Time step')

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
  %(prog)s exprs_midplane_s000000..000700.dat -x Psi_N -y ne --step 700
  %(prog)s exprs_midplane_s000000..000700.dat -x R -y currdens -o plot.png
  %(prog)s exprs_midplane_s000000..000700.dat --assume-si --n0 1.0
  %(prog)s exprs_midplane_s000000..000700.dat --jorek-units --n0 0.973 --mass 2.0

Available variables (depends on the expressions used in jorek2_postproc):
  R, ne, Psi_N, T_i, t, pres, currdens, Er, gradPdotCurv, dprof, zkiprof,
  ki_neo, mu_neo, T_e, eta_T, visco_T, J_bootstrap, JxB_R, gradP_R,
  vpar, V_ExB_pol, V_ExB_R, V_ExB_Z
        """)
    parser.add_argument('file', help='Path to midplane .dat file')
    parser.add_argument('-x', '--xvar', default='Psi_N',
                        help='Variable for x-axis (default: Psi_N)')
    parser.add_argument('-y', '--yvar', default='T_e',
                        help='Variable for y-axis (default: T_e)')
    parser.add_argument('--step', type=int, nargs='+',
                        help='Plot only specific time step(s)')
    parser.add_argument('--step-index', type=int, nargs='+',
                        help='Plot only specific step index/indices (0-based)')
    parser.add_argument('--all', action='store_true',
                        help='Plot all time steps (default: last step only)')
    parser.add_argument('--side', choices=['outer', 'inner', 'both'], default='outer',
                        help='Midplane side: outer (LFS, default), inner (HFS), or both')
    parser.add_argument('-o', '--output', help='Save plot to file')
    parser.add_argument('--title', help='Plot title')
    parser.add_argument('--no-show', action='store_true',
                        help='Do not display interactive window')
    parser.add_argument('--agg', action='store_true',
                        help='Use non-interactive Agg backend')

    # Unit conversion
    unit_group = parser.add_argument_group('Unit conversion')
    unit_group.add_argument('--assume-si', action='store_true', default=True,
                            help='Assume raw data is in SI units [default]')
    unit_group.add_argument('--jorek-units', action='store_true',
                            help='Raw data is in JOREK normalized units')
    unit_group.add_argument('--n0', type=float, default=0.973,
                            help='Central density (in 10^20 m^-3) [default: 0.973]')
    unit_group.add_argument('--mass', type=float, default=2.0,
                            help='Average ion mass in AMU [default: 2.0]')
    unit_group.add_argument('--cx', type=float,
                            help='Direct conversion factor for x-variable')
    unit_group.add_argument('--cy', type=float,
                            help='Direct conversion factor for y-variable')

    # List available variables
    parser.add_argument('--list-vars', action='store_true',
                        help='List available variables in the file and exit')

    args = parser.parse_args()

    if args.agg:
        matplotlib.use('Agg')

    # Parse the file
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
            print(f"  {v:20s} -> {info[3]:24s} [{unit}]")
        return

    assume_si = not args.jorek_units

    # Select steps: default = last step only
    selected_steps = None
    if args.step is not None:
        selected_steps = args.step
    elif args.step_index is not None:
        selected_steps = [time_steps[i]['step'] for i in args.step_index
                          if 0 <= i < len(time_steps)]
    elif not args.all:
        # Default: only the last time step
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
