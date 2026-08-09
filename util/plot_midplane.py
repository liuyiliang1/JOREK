#!/usr/bin/env python3
"""
Plot JOREK2 post-processed midplane data.

Read the exprs_midplane_*.dat (1D, toroidally averaged) or
exprs_midplane2d_*.dat (2D, full toroidal variation) files produced by
jorek2_postproc and plot any variable against any other, with appropriate
unit conversions.

Usage:
    python plot_midplane.py <midplane_file> [options]

Examples (1D — toroidally averaged profiles):
    # Plot T_e vs Psi_N for the last time step (default)
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat

    # Plot ne vs R for all time steps
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat -x R -y ne --all

    # Plot currdens vs Psi_N for a single time step
    python plot_midplane.py postproc/exprs_midplane_s000000..000700.dat -y currdens --step 500

Examples (2D — full toroidal variation, use with midplane2d output):
    # Overlay all toroidal angles (phi) as separate curves
    python plot_midplane.py postproc/exprs_midplane2d_*.dat -x R -y ne --2d

    # Plot a single toroidal plane (phi index 0)
    python plot_midplane.py postproc/exprs_midplane2d_*.dat -x R -y ne --2d --phi 0

    # 2D colormap: R vs phi, color = ne
    python plot_midplane.py postproc/exprs_midplane2d_*.dat -y ne --2d --pcolormesh

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
                # Fortran writes each name as character(len=23), fixed-width
                if not var_names:
                    line_raw = line.rstrip('\n\r')
                    body = line_raw[2:]  # skip '# '
                    var_names = []
                    for i in range(0, len(body), 23):
                        name = body[i:i+23].strip()
                        if name:
                            var_names.append(name)
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


def parse_midplane2d_file(filepath):
    """Parse a jorek2_postproc midplane2d output file.

    The midplane2d command evaluates expressions on the midplane at multiple
    toroidal angles.  The underlying ``res2d`` array has shape
    ``(nphi, npts, n_expr)``, meaning each ``# Column`` block corresponds to
    ONE toroidal angle and each row within a block contains the values of ALL
    expressions at one radial position.

    File format::

        # var1 var2 ... varN           (fixed 23-char fields)
        # time step #000100
          # Column 0001 of M           (M = nphi toroidal angles)
          v1_1  v1_2  ...  v1_N        (N = n_expr values, one per expression)
          v2_1  v2_2  ...  v2_N
          ...
          # Column 0002 of M
          ...

    Returns
    -------
    var_names : list of str
        Expression names from the header.
    time_steps : list of dict
        Each dict has keys 'step' (int), 'data' (dict: var_name -> 2D ndarray).
        Each 2D array has shape (npts, nphi).
    nphi : int
        Number of toroidal planes (deduced from Column count).
    """
    var_names = []
    time_steps = []
    current_step = None
    n_expr = 0                   # total expressions (from header line)
    nphi = 0                     # total toroidal angles (from Column markers)
    phi_idx = -1                 # current toroidal angle index (0-based)
    current_rows = []            # list of (n_expr,) rows for current phi slice
    all_phi_slices = []          # list of [rows_for_phi_0, rows_for_phi_1, ...]
    current_expr_name = None
    line_num = 0

    def _parse_fixed_width_header(line):
        """Parse the fixed-width header line into variable names.

        The Fortran code writes each expression name as character(len=23),
        so each name occupies exactly 23 characters in the header line.
        """
        line = line.rstrip('\n\r')
        if len(line) < 2:
            return []
        body = line[2:]   # skip '# '
        names = []
        for i in range(0, len(body), 23):
            name = body[i:i+23].strip()
            if name:
                names.append(name)
        return names

    def _flush_phi_slice():
        """Store current_rows as one toroidal slice."""
        nonlocal current_rows, phi_idx, all_phi_slices
        if phi_idx >= 0 and current_rows:
            # Pad all_phi_slices if needed
            while len(all_phi_slices) <= phi_idx:
                all_phi_slices.append([])
            all_phi_slices[phi_idx] = current_rows
        current_rows = []
        phi_idx = -1

    def _build_step_data():
        """Convert all_phi_slices into per-expression (npts, nphi) arrays."""
        if not all_phi_slices or n_expr == 0:
            return None
        # Determine npts from the first non-empty slice
        npts = 0
        for sl in all_phi_slices:
            if sl:
                npts = len(sl)
                break
        if npts == 0:
            return None
        nphi_actual = len(all_phi_slices)
        # Full 3D array: (nphi, npts, n_expr)
        full = np.full((nphi_actual, npts, n_expr), np.nan)
        for p in range(nphi_actual):
            rows = all_phi_slices[p]
            if len(rows) == npts:
                full[p, :, :] = np.array(rows)
            elif len(rows) > 0:
                # Partial slice – fill what we have
                full[p, :len(rows), :] = np.array(rows)
        # Reorganize to {varname: (npts, nphi)}
        data_dict = {}
        for e, vname in enumerate(var_names):
            if e < n_expr:
                data_dict[vname] = full[:, :, e].T  # shape (npts, nphi)
        # Also handle any placeholder names
        for e in range(len(var_names), n_expr):
            vname = f'__col_{e+1}__'
            data_dict[vname] = full[:, :, e].T
        return data_dict

    def _flush_step():
        """Save all accumulated toroidal slices as a completed time step."""
        nonlocal current_step, all_phi_slices, phi_idx, current_rows
        _flush_phi_slice()
        data_dict = _build_step_data()
        if current_step is not None and data_dict is not None:
            time_steps.append({
                'step': current_step,
                'data': data_dict
            })
        current_step = None
        all_phi_slices = []
        phi_idx = -1
        current_rows = []

    with open(filepath, 'r') as f:
        for line in f:
            line_num += 1
            stripped = line.strip()
            if not stripped:
                # Blank line signals end of time step
                _flush_step()
                continue

            if stripped.startswith('# time step'):
                # Flush previous step before starting a new one
                _flush_step()
                parts = stripped.split('#')
                try:
                    current_step = int(parts[-1].strip())
                except ValueError:
                    pass
                continue

            if stripped.startswith('# Column'):
                # Starting a new toroidal-angle block
                _flush_phi_slice()
                parts = stripped.split()
                try:
                    col_idx = int(parts[2])   # 1-based toroidal angle index
                    n_total = int(parts[4])   # total toroidal angles
                    phi_idx = col_idx - 1
                    nphi = n_total
                except (ValueError, IndexError):
                    phi_idx = -1
                continue

            if stripped.startswith('#'):
                # Header line with variable names ("# var1 var2 ...")
                if not var_names:
                    var_names = _parse_fixed_width_header(line)
                    n_expr = len(var_names)
                continue

            # Data line — n_expr values per row (all expressions at one radial
            # position, one toroidal angle)
            try:
                values = [float(x) for x in stripped.split()]
                if n_expr > 0 and len(values) == n_expr:
                    current_rows.append(values)
                elif n_expr == 0:
                    # Header not read yet; use this row to determine n_expr
                    current_rows.append(values)
            except ValueError:
                continue

    # Don't forget the last step
    _flush_step()

    # Determine nphi from data if header didn't provide it
    if nphi == 0 and time_steps:
        first_data = list(time_steps[0]['data'].values())[0]
        nphi = first_data.shape[1]

    return var_names, time_steps, nphi

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
# Midplane 2D plotting (toroidal variation preserved)
# ---------------------------------------------------------------------------

def plot_midplane2d(var_names, time_steps, nphi, x_var='R', y_var='T_e',
                    selected_steps=None, selected_phi=None, side='outer',
                    conversion_factor_x=None, conversion_factor_y=None,
                    assume_si=True, central_density=0.973, central_mass=2.0,
                    as_pcolormesh=False,
                    title=None, output_file=None, interactive=True):
    """Plot 2D midplane data with full toroidal variation.

    Parameters
    ----------
    var_names : list of str
        Expression names.
    time_steps : list of dict
        Parsed 2D time step data (data[var_name] is 2D: npts x nphi).
    nphi : int
        Number of toroidal planes.
    x_var, y_var : str
        Variable names for x and y axes.
    selected_steps : list of int or None
        Which time steps to plot. None = last step only.
    selected_phi : int, list of int, or None
        Which toroidal plane(s) to plot. None = all phi.
    conversion_factor_x, conversion_factor_y : float or None
        Direct scaling factors.
    assume_si : bool
        If True, assume raw data is in SI.
    as_pcolormesh : bool
        If True, use pcolormesh (2D colormap) instead of overlay curves.
    title : str or None
        Plot title.
    output_file : str or None
        Save to this file instead of showing.
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

    # Select time steps
    if selected_steps is not None and len(selected_steps) > 0:
        steps_to_plot = [ts for ts in time_steps if ts['step'] in selected_steps]
    else:
        steps_to_plot = list(time_steps)   # None or empty = all steps

    if not steps_to_plot:
        print("ERROR: No matching time steps found.")
        sys.exit(1)

    # ---- Side filtering (outer / inner midplane) ----
    if side in ('outer', 'inner') and 'Psi_N' in var_names:
        for ts in steps_to_plot:
            psi = ts['data']['Psi_N']   # (npts, nphi)
            # Use first phi column to find magnetic axis (min Psi_N)
            psi_col0 = psi[:, 0]
            idx_min = np.argmin(psi_col0)
            if side == 'outer':
                # Low-field side: from magnetic axis outward
                ts['data'] = {k: v[idx_min:, :] for k, v in ts['data'].items()}
            else:
                # High-field side: from inner wall to magnetic axis
                ts['data'] = {k: v[:idx_min + 1, :] for k, v in ts['data'].items()}

    # Determine which phi indices to plot
    if selected_phi is not None:
        if isinstance(selected_phi, int):
            phi_indices = [selected_phi]
        else:
            phi_indices = list(selected_phi)
    else:
        phi_indices = list(range(nphi))

    # Validate phi indices
    phi_indices = [p for p in phi_indices if 0 <= p < nphi]
    if not phi_indices:
        print(f"ERROR: No valid phi indices (nphi={nphi}).")
        sys.exit(1)

    # Build phi labels (in degrees)
    phi_deg = [360.0 * p / nphi for p in phi_indices]

    # ---- Pcolormesh mode (2D colormap: φ vs Psi_N) ----
    if as_pcolormesh:
        n_steps = len(steps_to_plot)
        ncols = min(n_steps, 3)
        nrows = (n_steps + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(6.3 * ncols, 5.0 * nrows),
                                 squeeze=False)
        full_phi_deg = np.array([360.0 * p / nphi for p in range(nphi)])

        for idx, ts in enumerate(steps_to_plot):
            ax = axes[idx // ncols][idx % ncols]
            y_data = ts['data'][y_var]     # (npts, nphi)

            y_disp, y_unit = convert_variable(
                y_data, y_var, conversion_factor_y,
                central_density, central_mass, assume_si)

            # Use Psi_N as y-axis (flux coordinate), or radial index if unavailable
            if 'Psi_N' in ts['data']:
                psi_2d = ts['data']['Psi_N']                    # (npts, nphi)
                # Psi_N should be nearly constant along phi; use mean
                psi_y = np.mean(psi_2d, axis=1)                 # (npts,)
                y_label = r'$\Psi_N$'
            else:
                psi_y = np.arange(y_data.shape[0])
                y_label = 'Radial index'

            # pcolormesh: x=phi, y=Psi_N, color=y_var
            im = ax.pcolormesh(full_phi_deg, psi_y, y_disp,
                               shading='auto', cmap='viridis')
            plt.colorbar(im, ax=ax, label=_build_axis_label(y_var, y_unit))
            ax.set_xlabel(r'$\phi$ (deg)')
            ax.set_ylabel(y_label)
            ax.set_title(f"step {ts['step']}")
        # Hide unused subplots
        for idx in range(n_steps, nrows * ncols):
            axes[idx // ncols][idx % ncols].set_visible(False)
        fig.tight_layout()

    # ---- Overlay mode (curves for each phi) ----
    else:
        fig, ax = plt.subplots(figsize=(12.6, 10.0))
        fig.subplots_adjust(left=0.12, bottom=0.12, right=0.94, top=0.94)
        ax.set_box_aspect(1.15 / 1.45)

        n_steps = len(steps_to_plot)
        n_phis = len(phi_indices)
        n_curves = n_phis * n_steps

        # Color maps
        time_cmap = plt.cm.plasma       # distinct colors for time steps
        phi_cmap = plt.cm.viridis       # for phi colorbar fallback

        # Line style cycling for phi (when both dimensions vary)
        phi_linestyles = ['-', '--', '-.', ':']
        phi_alphas = [1.0, 0.75, 0.55, 0.4]

        for si, ts in enumerate(steps_to_plot):
            x_data = ts['data'][x_var]     # (npts, nphi)
            y_data = ts['data'][y_var]     # (npts, nphi)

            x_disp, x_unit = convert_variable(
                x_data, x_var, conversion_factor_x, central_density, central_mass, assume_si)
            y_disp, y_unit = convert_variable(
                y_data, y_var, conversion_factor_y, central_density, central_mass, assume_si)

            # Base color from time step position
            t_color = time_cmap(float(si) / max(n_steps - 1, 1))

            for pi, p in enumerate(phi_indices):
                x_slice = x_disp[:, p]
                y_slice = y_disp[:, p]

                if n_phis > 1 and n_steps > 1:
                    # Both dimensions vary: color by time, vary style by phi
                    color = t_color
                    ls = phi_linestyles[pi % len(phi_linestyles)]
                    alpha = phi_alphas[min(pi, len(phi_alphas) - 1)]
                    label = f"step {ts['step']}, φ={phi_deg[pi]:.0f}°"
                elif n_phis > 1:
                    # Single time step, multiple phi: color by phi
                    color = phi_cmap(float(p) / max(nphi - 1, 1))
                    ls = '-'
                    alpha = 1.0
                    label = f"φ={phi_deg[pi]:.0f}°"
                else:
                    # Single phi, multiple time steps: color by time
                    color = t_color
                    ls = '-'
                    alpha = 1.0
                    label = f"step {ts['step']}"

                ax.plot(x_slice, y_slice, color=color, linestyle=ls,
                        alpha=alpha, label=label, linewidth=2)

        ax.set_xlabel(_build_axis_label(x_var, x_unit))
        ax.set_ylabel(_build_axis_label(y_var, y_unit))

        if n_curves <= 15:
            ax.legend(loc='best', fontsize=14)
        elif n_phis > 1 and n_steps > 1:
            # Too many to legend; add time-step colorbar
            sm = plt.cm.ScalarMappable(cmap=time_cmap,
                                       norm=plt.Normalize(
                                           vmin=steps_to_plot[0]['step'],
                                           vmax=steps_to_plot[-1]['step']))
            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_label('Time step')
        else:
            sm = plt.cm.ScalarMappable(cmap=phi_cmap,
                                       norm=plt.Normalize(vmin=0, vmax=360))
            cbar = plt.colorbar(sm, ax=ax)
            cbar.set_label(r'$\phi$ (deg)')

        ax.grid(True, alpha=0.3)

    if title:
        fig.suptitle(title) if as_pcolormesh else ax.set_title(title)

    if output_file:
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved to {output_file}")
    if interactive:
        plt.show()


def _build_axis_label(var_name, unit_label):
    """Build a LaTeX axis label from variable info and unit string."""
    info = VARIABLE_INFO.get(var_name, (var_name, '', 1.0, var_name))
    display_name = info[3]
    if unit_label:
        return f'{display_name} ({unit_label})'
    return display_name


# ---------------------------------------------------------------------------
# Plotting (1D)
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
    parser.add_argument('--2d', dest='is_2d', action='store_true',
                        help='Parse file as midplane2d output (R-phi matrix, full toroidal variation)')
    parser.add_argument('--phi', type=int, nargs='+',
                        help='Toroidal plane index(es) to plot, 0-based '
                             '(only for --2d overlay).  Use --phi-deg for degrees.')
    parser.add_argument('--phi-deg', type=float, nargs='+', dest='phi_deg_input',
                        help='Toroidal angle(s) in degrees (only for --2d overlay). '
                             'E.g. --phi-deg 0 90 180 270')
    parser.add_argument('--pcolormesh', action='store_true',
                        help='Use 2D colormap (R vs phi) instead of overlay curves (only for --2d)')
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

    # ---- 2D mode (midplane2d) ----
    if args.is_2d:
        var_names, time_steps, nphi = parse_midplane2d_file(args.file)

        if not var_names:
            print("ERROR: Could not parse variable names from header.")
            print("First 5 lines of file:")
            with open(args.file) as f:
                for i, l in enumerate(f):
                    if i >= 5:
                        break
                    print(f"  {i+1}: {l.rstrip()}")
            sys.exit(1)
        if not time_steps:
            print("ERROR: No time step data found in file.")
            print(f"Header parsed: {len(var_names)} variables: {var_names}")
            print("First 10 lines of file:")
            with open(args.file) as f:
                for i, l in enumerate(f):
                    if i >= 10:
                        break
                    print(f"  {i+1}: {l.rstrip()}")
            sys.exit(1)

        first_data = list(time_steps[0]['data'].values())[0]
        npts, nphi_actual = first_data.shape
        print(f"File: {args.file}")
        print(f"Variables: {var_names}")
        print(f"Time steps: {len(time_steps)} steps "
              f"({time_steps[0]['step']} to {time_steps[-1]['step']})")
        print(f"Grid: {npts} radial × {nphi_actual} toroidal points")

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
            selected_steps = [time_steps[-1]['step']]
            print(f"Plotting last step only (use --all for all {len(time_steps)} steps)")

        # Phi selection
        if args.phi_deg_input is not None:
            # Convert degrees to 0-based indices
            selected_phi = [int(round(d * nphi / 360.0)) % nphi
                            for d in args.phi_deg_input]
        else:
            selected_phi = args.phi  # None = all phi, or list of indices

        plot_midplane2d(
            var_names, time_steps, nphi,
            x_var=args.xvar, y_var=args.yvar,
            selected_steps=selected_steps,
            selected_phi=selected_phi,
            side=args.side,
            conversion_factor_x=args.cx,
            conversion_factor_y=args.cy,
            assume_si=assume_si,
            central_density=args.n0,
            central_mass=args.mass,
            as_pcolormesh=args.pcolormesh,
            title=args.title,
            output_file=args.output,
            interactive=not args.no_show,
        )
        return

    # ---- 1D mode (midplane, toroidally averaged) ----
    var_names, time_steps = parse_midplane_file(args.file)

    if not var_names:
        print("ERROR: Could not parse variable names from header.")
        print("First 5 lines of file:")
        with open(args.file) as f:
            for i, l in enumerate(f):
                if i >= 5:
                    break
                print(f"  {i+1}: {l.rstrip()}")
        sys.exit(1)
    if not time_steps:
        print("ERROR: No time step data found in file.")
        print(f"Header parsed: {len(var_names)} variables: {var_names}")
        print("First 10 lines of file:")
        with open(args.file) as f:
            for i, l in enumerate(f):
                if i >= 10:
                    break
                print(f"  {i+1}: {l.rstrip()}")
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
