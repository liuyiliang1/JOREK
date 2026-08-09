#!/usr/bin/env python3
"""
Generate stpts file for jorek2_poincare and plot the resulting Poincaré plot.

Usage:
  python generate_poincare.py              # Generate stpts only
  python generate_poincare.py --plot       # Generate stpts and plot poinc_R-Z.dat
  python generate_poincare.py --plot-only  # Only plot poinc_R-Z.dat (skip stpts generation)

Workflow:
  1. python generate_poincare.py           # Create stpts file
  2. ./jorek2_poincare                     # Run field line tracing
  3. python generate_poincare.py --plot    # Plot the results
"""

import numpy as np
import sys
import os

# ==============================================================================
# Equilibrium data (extracted from logfile / special_equilibrium_points.dat)
# ==============================================================================
R_axis  = 1.24356       # Magnetic axis R
Z_axis  = -0.00496      # Magnetic axis Z
Psi_axis = 0.59915      # Psi at axis
Psi_bnd  = 0.12573      # Psi at separatrix (X-point)

# X-points
R_xp_lower = 0.85724
Z_xp_lower = -1.34359
R_xp_upper = 0.85726
Z_xp_upper = 1.34355

# LCFS shape parameters (from logfile)
R_geo    = 1.09218      # Geometric major radius
Z_geo    = -0.00021     # Geometric Z center
a_min    = 0.59562      # Minor radius
kappa    = 2.23451      # Elongation
delta_U  = 0.38282      # Upper triangularity
delta_L  = 0.38150      # Lower triangularity

# Inner midplane LCFS point
R_inner_mid = 0.49127

# ==============================================================================
# Estimate outboard midplane separatrix position
# For LCFS: R = R_geo + a*cos(theta + delta*sin(theta))
# At outboard midplane (theta=0): R_out = R_geo + a_min
# ==============================================================================
R_outer_mid = R_geo + a_min  # ≈ 1.6878

def print_info():
    """Print equilibrium summary."""
    print("=" * 60)
    print("  Poincare stpts Generator")
    print("=" * 60)
    print(f"  Magnetic axis:    R={R_axis:.5f}, Z={Z_axis:.5f}")
    print(f"  Psi_axis:         {Psi_axis:.5f}")
    print(f"  Psi_separatrix:   {Psi_bnd:.5f}")
    print(f"  Inner midplane:   R={R_inner_mid:.5f}")
    print(f"  Outer midplane:   R={R_outer_mid:.5f} (estimated)")
    print(f"  X-points:         R={R_xp_lower:.5f}, Z=±{abs(Z_xp_lower):.5f}")
    print(f"  LCFS: κ={kappa:.3f}, δ_U={delta_U:.3f}, a={a_min:.5f}")
    print("=" * 60)


# ==============================================================================
# Generate stpts file
# ==============================================================================
def generate_stpts(n_lines=60, n_turns=500, filename="stpts"):
    """
    Generate stpts file for jorek2_poincare.

    Starting points are placed along the outboard midplane (Z = Z_axis),
    from near the magnetic axis to near the separatrix.

    Parameters
    ----------
    n_lines : int
        Number of field lines to trace.
    n_turns : int
        Number of poloidal turns per field line.
    filename : str
        Output filename.
    """
    # Use non-uniform spacing: more points near separatrix (interesting region)
    # psi_norm from ~0.05 (near axis) to ~0.98 (near separatrix)
    psi_norm_start = np.linspace(0.0, 1.0, n_lines)**0.5 * 0.98 + 0.02

    # Approximate mapping: psi_n → R at outboard midplane
    # Use a simple model: R ≈ R_axis + (R_outer_mid - R_axis) * psi_n^(something)
    # Actually, for a typical equilibrium, sqrt(psin) ~ r/a, and R ~ R_axis + r
    # Let's use: R_start = R_axis + (R_outer_mid - R_axis) * sqrt(psi_n)
    R_start = R_axis + (R_outer_mid - R_axis) * np.sqrt(psi_norm_start)

    Z_start = np.full(n_lines, Z_axis)  # All at midplane Z
    phi_start = np.zeros(n_lines)        # Toroidal angle 0

    with open(filename, 'w') as f:
        f.write("# Poincare plot starting points (generated automatically)\n")
        f.write(f" {n_lines}\n")
        f.write("# nr   R_start        Z_start        phi_start   n_turns\n")

        # Write grouped points: first point at index 1, last at n_lines
        # Points in between are linearly interpolated by jorek2_poincare
        f.write(f" {1:4d}  {R_start[0]:.6f}  {Z_start[0]:.6f}  {phi_start[0]:.6f}  {n_turns}\n")
        f.write(f" {n_lines:4d}  {R_start[-1]:.6f}  {Z_start[-1]:.6f}  {phi_start[-1]:.6f}  {n_turns}\n")

    print(f"\n  Generated: {filename}")
    print(f"    {n_lines} field lines, {n_turns} turns each")
    print(f"    R range: {R_start[0]:.5f} → {R_start[-1]:.5f}")
    print(f"    Z start: {Z_start[0]:.5f}")

    return R_start, Z_start


def generate_stpts_fine(filename="stpts"):
    """
    Alternative: Generate stpts with explicit points grouped into blocks
    for more control over spacing. Uses the nr interpolation feature.

    This creates denser spacing near the separatrix and axis.
    """
    # Strategy: create multiple groups, each with linear interpolation
    # Group 1: Near axis (nr=1 to nr=15)
    # Group 2: Middle region (nr=15 to nr=45)
    # Group 3: Near separatrix (nr=45 to nr=60)

    n_total = 200
    n_turns = 2500

    # Fine spacing near separatrix, coarser in core
    psi_breakpoints = [
        (1,   0.00),   # nr, psi_n
        (200,  0.99),
    ]

    lines = []
    lines.append("# Poincare plot starting points (fine spacing near separatrix)")
    lines.append(f" {n_total}")
    lines.append("# nr   R_start        Z_start        phi_start   n_turns")

    for nr, psi_n in psi_breakpoints:
        R = R_axis + (R_outer_mid - R_axis) * np.sqrt(psi_n)
        lines.append(f" {nr:4d}  {R:.6f}  {Z_axis:.6f}  0.000000  {n_turns}")

    with open(filename, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    print(f"\n  Generated (fine): {filename}")
    print(f"    {n_total} field lines, {n_turns} turns each")
    return


# ==============================================================================
# Plotting
# ==============================================================================
def plot_poincare_psi_theta(datafile="poinc_rho-theta.dat", output="poincare_psi_theta.png",
                            q_surfaces=None):
    """
    Plot the Poincaré plot in flux coordinates (psi_n vs poloidal angle theta).

    Each field line is assigned a distinct color.
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    plt.rcParams.update({
        'font.size': 28,
        'font.weight': 'bold',
        'axes.linewidth': 3,
    })

    field_lines = read_poincare_psi_theta(datafile)

    if field_lines is None or len(field_lines) == 0:
        print(f"ERROR: No data found in {datafile}")
        print("Run ./jorek2_poincare first to generate the data.")
        return

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Colormap per field line
    n_lines = len(field_lines)
    cmap = mpl.cm.get_cmap('tab10')
    colors = cmap(np.linspace(0, 1, min(n_lines, 10)))

    psi_all = []
    for i, (rho, theta) in enumerate(field_lines):
        psi_n = rho**2
        theta_norm = np.mod(theta, 2 * np.pi)
        ax.scatter(psi_n, theta_norm, s=0.5, color=colors[i % len(colors)],
                   alpha=0.6, rasterized=True)
        psi_all.append(psi_n)

    psi_all = np.concatenate(psi_all)

    ax.set_xlabel(r'$\psi_{nor}$')
    ax.set_ylabel(r'$\theta$ [rad]')
    # Adaptive x-axis
    psi_min = max(np.nanmin(psi_all), 0.0)
    psi_max = np.nanmax(psi_all)
    psi_margin = 0.05 * (psi_max - psi_min)
    ax.set_xlim(psi_min - psi_margin, psi_max + psi_margin)
    ax.set_ylim(0, 2 * np.pi)
    ax.tick_params(width=3, length=8)
    ax.grid(True, alpha=0.3)

    # Draw rational surface lines
    if q_surfaces is None:
        q_surfaces = {2: 0.56, 3: 0.67, 4: 0.74}
    for q_val, psi_n_q in q_surfaces.items():
        if 0 < psi_n_q < 1:
            ax.axvline(x=psi_n_q, color='black',
                      linestyle='--', alpha=0.7, linewidth=4)
            ax.text(psi_n_q, -0.22, f'q={q_val}', color='black',
                    fontsize=24, fontweight='bold', ha='center', va='top',
                    transform=ax.get_xaxis_transform())

    plt.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches='tight')
    print(f"\n  Psi-theta plot saved to: {output}")
    plt.show()


def read_poincare_psi_theta(filename="poinc_rho-theta.dat"):
    """
    Read jorek2_poincare psi_n-theta output file.

    Returns a list of (rho, theta) arrays, one per field line.
    Field lines are separated by blank lines in the data file.
    """
    field_lines = []
    rho_cur, theta_cur = [], []

    if not os.path.exists(filename):
        print(f"WARNING: {filename} not found.")
        return None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                if rho_cur:
                    field_lines.append((np.array(rho_cur), np.array(theta_cur)))
                    rho_cur, theta_cur = [], []
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    rho_cur.append(float(parts[0]))
                    theta_cur.append(float(parts[1]))
                except ValueError:
                    continue

    if rho_cur:
        field_lines.append((np.array(rho_cur), np.array(theta_cur)))

    return field_lines


def plot_poincare(datafile="poinc_R-Z.dat", output="poincare_plot.png"):
    """
    Plot the Poincaré plot from jorek2_poincare output.

    Each field line is assigned a distinct color from a colormap,
    similar to plot_poincare.py.
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    plt.rcParams.update({
        'font.size': 28,
        'font.weight': 'bold',
        'axes.linewidth': 3,
    })

    # Read equilibrium points
    eq_points = read_equilibrium_points()

    # Read poincare data: list of (R, Z) per field line
    field_lines = read_poincare_RZ(datafile)

    if field_lines is None or len(field_lines) == 0:
        print(f"ERROR: No data found in {datafile}")
        print("Run ./jorek2_poincare first to generate the data.")
        return

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(8, 12))
    ax.set_aspect('equal')

    # Colormap: cycle through colors per field line
    n_lines = len(field_lines)
    cmap = mpl.cm.get_cmap('tab10')
    colors = cmap(np.linspace(0, 1, min(n_lines, 10)))

    for i, (R, Z) in enumerate(field_lines):
        ax.scatter(R, Z, s=0.5, color=colors[i % len(colors)],
                   alpha=0.6, rasterized=True)

    # Mark magnetic axis
    ax.plot(eq_points['R_axis'], eq_points['Z_axis'], 'r+',
            markersize=15, markeredgewidth=2, label='Magnetic Axis')

    # Mark X-points
    ax.plot(eq_points['R_xp_lower'], eq_points['Z_xp_lower'], 'bx',
            markersize=12, markeredgewidth=2, label='Lower X-point')
    ax.plot(eq_points['R_xp_upper'], eq_points['Z_xp_upper'], 'bx',
            markersize=12, markeredgewidth=2, label='Upper X-point')

    # Mark boundary point
    ax.plot(eq_points['R_bnd'], eq_points['Z_bnd'], 'g^',
            markersize=8, label='Boundary (LCFS)')

    # Mark strike points if available
    for i, (r, z) in enumerate(eq_points['strike_points']):
        label = 'Strike points' if i == 0 else None
        ax.plot(r, z, 'ms', markersize=6, label=label)

    ax.set_xlabel('R [m]')
    ax.set_ylabel('Z [m]')
    ax.tick_params(width=3, length=8)
    ax.legend(loc='upper right', fontsize=14, framealpha=0.8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(output, dpi=150, bbox_inches='tight')
    print(f"\n  Plot saved to: {output}")
    plt.show()


def read_equilibrium_points(filename="special_equilibrium_points.dat"):
    """Parse special_equilibrium_points.dat (ASCII format)."""
    points = {
        'R_axis': R_axis,
        'Z_axis': Z_axis,
        'R_xp_lower': R_xp_lower,
        'Z_xp_lower': Z_xp_lower,
        'R_xp_upper': R_xp_upper,
        'Z_xp_upper': Z_xp_upper,
        'R_bnd': R_xp_lower,   # Boundary = lower X-point for this case
        'Z_bnd': Z_xp_lower,
        'strike_points': [],
    }

    if os.path.exists(filename):
        with open(filename, 'r') as f:
            content = f.read()

        lines = content.split('\n')
        current_section = None
        for line in lines:
            line = line.strip()
            if line.startswith('# Magnetic Axis'):
                current_section = 'axis'
            elif line.startswith('# Lower X-Point'):
                current_section = 'xp_lower'
            elif line.startswith('# Upper X-Point'):
                current_section = 'xp_upper'
            elif line.startswith('# Boundary'):
                current_section = 'bnd'
            elif line.startswith('# Strike point'):
                current_section = 'strike'
            elif line.startswith('# Inner midplane'):
                current_section = 'inner_mid'
            elif line.startswith('# Limiter'):
                current_section = 'limiter'
            elif line and not line.startswith('#'):
                # Parse numbers
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        r, z = float(parts[0]), float(parts[1])
                        if current_section == 'axis':
                            points['R_axis'] = r
                            points['Z_axis'] = z
                        elif current_section == 'xp_lower':
                            points['R_xp_lower'] = r
                            points['Z_xp_lower'] = z
                        elif current_section == 'xp_upper':
                            points['R_xp_upper'] = r
                            points['Z_xp_upper'] = z
                        elif current_section == 'bnd':
                            points['R_bnd'] = r
                            points['Z_bnd'] = z
                        elif current_section == 'strike':
                            points['strike_points'].append((r, z))
                    except ValueError:
                        pass

    return points


def read_poincare_RZ(filename="poinc_R-Z.dat"):
    """
    Read jorek2_poincare output file.

    Returns a list of (R, Z) arrays, one per field line.
    Field lines are separated by blank lines in the data file.
    """
    field_lines = []
    R_cur, Z_cur = [], []

    if not os.path.exists(filename):
        print(f"WARNING: {filename} not found.")
        return None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                if R_cur:
                    field_lines.append((np.array(R_cur), np.array(Z_cur)))
                    R_cur, Z_cur = [], []
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    R_cur.append(float(parts[0]))
                    Z_cur.append(float(parts[1]))
                except ValueError:
                    continue

    if R_cur:
        field_lines.append((np.array(R_cur), np.array(Z_cur)))

    return field_lines


# ==============================================================================
# Main
# ==============================================================================
if __name__ == "__main__":
    plot_only = '--plot-only' in sys.argv
    do_plot = '--plot' in sys.argv or plot_only

    # Determine plot type: rz, psi, or both (default)
    plot_type = 'both'
    for i, arg in enumerate(sys.argv):
        if arg in ('--plot', '--plot-only'):
            if i + 1 < len(sys.argv) and sys.argv[i+1] in ('rz', 'psi', 'both'):
                plot_type = sys.argv[i+1]

    print_info()

    if not plot_only:
        generate_stpts_fine("stpts")
        print("\n" + "=" * 60)
        print("  Next steps:")
        print("    1. Run:  ./jorek2_poincare")
        print("    2. Plot: python generate_poincare.py --plot [rz|psi|both]")
        print("=" * 60)

    if do_plot:
        if plot_type in ('rz', 'both'):
            plot_poincare()
        if plot_type in ('psi', 'both'):
            plot_poincare_psi_theta()
