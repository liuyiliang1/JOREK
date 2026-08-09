#!/usr/bin/env python3
"""
Plot target plate profiles from target_strike_profile.dat output of jorek2_target2vtkn.

Usage:
    python plot_target_strike.py <data_dir> [options]

Output columns in target_strike_profile.dat:
    plate_id  strike_dist  R  Z  angle  density  Ti  Te  Vpar  nv  nvT_gam  KparT  psi  sputter
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os
import argparse

# =============================================================================
# Sputtering yield calculation (ported from mod_eckstein_thompson.f90)
# =============================================================================

def eckstein_sputter_yield(Z_ion, Z_target, M_ion, M_target, E_eV, theta_deg):
    """
    Physical sputtering yield using Eckstein formula.

    Parameters
    ----------
    Z_ion, Z_target : int
        Atomic numbers of projectile and target
    M_ion, M_target : float
        Masses (amu) of projectile and target
    E_eV : float or array
        Incident energy in eV
    theta_deg : float or array
        Incident angle in degrees (0 = normal incidence)

    Returns
    -------
    yield_val : float or array
        Sputtering yield (atoms/ion)
    """
    # Surface binding energy
    Es_dict = {74: 8.9931, 6: 7.42, 42: 6.83, 26: 4.29, 28: 4.44, 29: 3.49, 13: 3.36, 14: 4.70}
    Es = Es_dict.get(Z_target, 0.5 * M_target / 1000.0)

    # Threshold energy
    if Z_ion == 1:
        Eth_dict = {74: 220.5, 6: 27.64, 42: 180.0, 26: 120.0}
        Eth = Eth_dict.get(Z_target, 4.0 * Es * (M_ion + M_target)**2 / (4.0 * M_ion * M_target))
    elif Z_ion == 2:
        Eth_dict = {74: 62.06, 6: 52.98}
        Eth = Eth_dict.get(Z_target, 4.0 * Es * (M_ion + M_target)**2 / (4.0 * M_ion * M_target))
    else:
        Eth = 4.0 * Es * (M_ion + M_target)**2 / (4.0 * M_ion * M_target)

    # Surface coordination number
    ns_dict = {74: 0.06325, 6: 0.11286, 42: 0.075, 26: 0.085}
    ns = ns_dict.get(Z_target, 0.1)

    # Lindhard screening length
    aL = 0.4685 / np.sqrt(Z_ion**(2/3) + Z_target**(2/3))

    # Thomas-Fermi energy
    ETF = Z_ion * Z_target * 14.4 / aL * (M_ion + M_target) / M_target

    # Energy transfer factor
    gamma1 = 4.0 * M_ion * M_target / (M_ion + M_target)**2

    # Sputtering efficiency Q
    Q_val = (1.633 * Z_ion**(2/3) * Z_target**(2/3) *
             (Z_ion**(2/3) + Z_target**(2/3))**(1/3) *
             M_ion**(5/6) * M_target**(1/6) / (M_ion + M_target) *
             (0.15 + 0.05 * (M_target / M_ion)) /
             (1.0 + 0.05 * (M_target / M_ion)**1.6) * Es**(-2/3))

    # Surface roughness parameter f
    f_val = np.sqrt(Es) * (0.94 - 1.33e-3 * M_target / M_ion)

    # Normalized energy
    epsilon1 = np.maximum(E_eV, 1e-10) / ETF

    # Nuclear stopping power (KrC)
    w = epsilon1 + 0.1728 * np.sqrt(epsilon1) + 0.008 * epsilon1**0.1504
    Sn_KRC = 0.5 * np.log(1.0 + 1.2288 * epsilon1) / w

    # Normal incidence yield
    ratio = np.maximum(Eth / E_eV, 0.0)
    mask = (E_eV > Eth)
    Yn = np.where(mask, Q_val * Sn_KRC * (1.0 - ratio**(2/3)) * (1.0 - ratio)**2, 0.0)

    # Angle-dependent theta_star
    aopt = (np.pi/2 - aL * ns**(1/3) *
            (2.0 * epsilon1 * np.sqrt(Es / (gamma1 * np.maximum(E_eV, 1e-10))))**(-0.5))
    theta_star_rad = np.where(mask, aopt, np.pi/2)

    # Angle dependence
    theta_rad = np.abs(theta_deg) * np.pi / 180.0
    cos_theta = np.maximum(np.cos(theta_rad), 0.01)

    angle_factor = (cos_theta**(-f_val) *
                    np.exp(f_val * (1.0 - 1.0/cos_theta) *
                           np.cos(theta_star_rad)))

    return Yn * angle_factor


def chemical_sputtering(T_target_eV, E_in_eV, flux_in):
    """
    Chemical sputtering yield (H on C only, returns 0 for other targets).
    Ported from mod_eckstein_thompson.f90.

    Parameters
    ----------
    T_target_eV : float or array
        Target surface temperature in eV
    E_in_eV : float or array
        Incident ion energy in eV
    flux_in : float or array
        Incident particle flux in m^-2 s^-1
    """
    E_TF_H = 415.0
    Q_H = 0.035
    D_H = 250.0
    E_dam = 15.0
    E_des = 2.0
    E_rel = 1.8
    E_therm = 1.7

    ratio_vals = np.maximum(E_in_eV / E_TF_H, 1e-10)
    Sn_E0 = (0.5 * np.log(1.0 + 1.2288 * ratio_vals) /
             (ratio_vals + 0.1728 * np.sqrt(ratio_vals) + 0.008 * ratio_vals**0.1504))

    Ydam = Q_H * Sn_E0 * (1.0 - (E_dam / np.maximum(E_in_eV, E_dam+1))**(2/3)) * (1.0 - E_dam / np.maximum(E_in_eV, E_dam+1))**2
    Ydes = Q_H * Sn_E0 * (1.0 - (E_des / np.maximum(E_in_eV, E_des+1))**(2/3)) * (1.0 - E_des / np.maximum(E_in_eV, E_des+1))**2

    T_K = T_target_eV * 11604.5  # eV to K
    flux_judge = 1e30 * np.exp(-1.4 / np.maximum(T_target_eV, 0.01))

    c = 1.0 / (1.0 + 3e-23 * np.where(flux_in > flux_judge, flux_in, 0.0) +
               3e7 * np.exp(-1.4 / np.maximum(T_target_eV, 0.01)) * np.where(flux_in <= flux_judge, 1.0, 0.0))

    c_sp3 = (c * (2e-32 * flux_in + np.exp(-E_therm / np.maximum(T_target_eV, 0.01))) /
             (2e-32 * flux_in + (1.0 + 2e29 / np.maximum(flux_in, 1.0) * np.exp(-E_rel / np.maximum(T_target_eV, 0.01))) *
              np.exp(-E_therm / np.maximum(T_target_eV, 0.01))))

    Ytherm = (c_sp3 * 0.033 * np.exp(-E_therm / np.maximum(T_target_eV, 0.01)) /
              (2e-32 * flux_in + np.exp(-E_therm / np.maximum(T_target_eV, 0.01))))
    Ysurf = c_sp3 * Ydes / (1.0 + np.exp((E_in_eV - 65.0) / 40.0))

    Ydam = np.maximum(Ydam, 0.0)
    Ydes = np.maximum(Ydes, 0.0)
    Ytherm = np.maximum(Ytherm, 0.0)
    Ysurf = np.maximum(Ysurf, 0.0)

    return Ytherm * (1.0 + D_H * Ydam) + Ysurf


# =============================================================================
# Data reading
# =============================================================================

def read_target_data(filepath):
    """Read target_strike_profile.dat and return dict of plates."""
    plates = {}
    current_plate = None
    current_data = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                if current_plate is not None and current_data:
                    plates[current_plate] = np.array(current_data)
                    current_data = []
                    current_plate = None
                continue
            if line.startswith('# Plate'):
                # Save previous plate if any
                if current_plate is not None and current_data:
                    plates[current_plate] = np.array(current_data)
                    current_data = []
                # Parse plate id
                parts = line.split()
                current_plate = int(parts[2])
            elif line.startswith('#'):
                continue
            else:
                vals = [float(x) for x in line.split()]
                current_data.append(vals)

        # Don't forget the last plate
        if current_plate is not None and current_data:
            plates[current_plate] = np.array(current_data)

    # Sort each plate by strike_dist and remove duplicates
    for pid in plates:
        data = plates[pid]
        # Remove duplicate rows (same strike_dist, R, Z)
        _, unique_idx = np.unique(np.round(data[:, 1:4], decimals=6), axis=0, return_index=True)
        data = data[np.sort(unique_idx)]
        # Sort by strike_dist for monotonic x-axis
        data = data[data[:, 1].argsort()]
        plates[pid] = data

    return plates


def get_xlim(sdist_cm):
    """Adaptive x-axis: 1.2x of SOL and private region ranges."""
    pos = sdist_cm[sdist_cm >= 0]
    neg = sdist_cm[sdist_cm < 0]
    xmax = np.max(pos) * 1.2 if len(pos) > 0 else 1.0
    xmin = np.min(neg) * 1.2 if len(neg) > 0 else -1.0
    xmax = max(xmax, 1.0)   # at least 1 cm
    xmin = min(xmin, -0.5)  # at least -0.5 cm
    return xmin, xmax


# =============================================================================
# Plotting
# =============================================================================

def setup_style():
    """Configure matplotlib to match MATLAB figure style."""
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
    rcParams['font.size'] = 24
    rcParams['font.weight'] = 'bold'
    rcParams['axes.linewidth'] = 4
    rcParams['lines.linewidth'] = 4
    rcParams['lines.solid_capstyle'] = 'round'
    rcParams['xtick.major.width'] = 2
    rcParams['ytick.major.width'] = 2
    rcParams['xtick.major.size'] = 8
    rcParams['ytick.major.size'] = 8
    rcParams['xtick.minor.size'] = 4
    rcParams['ytick.minor.size'] = 4
    rcParams['xtick.labelsize'] = 28
    rcParams['ytick.labelsize'] = 28
    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.direction'] = 'in'
    rcParams['legend.fontsize'] = 24
    rcParams['legend.frameon'] = False
    rcParams['axes.labelsize'] = 32
    rcParams['axes.labelweight'] = 'bold'
    rcParams['axes.titlesize'] = 28
    rcParams['figure.facecolor'] = 'white'
    rcParams['axes.facecolor'] = 'white'
    rcParams['savefig.dpi'] = 300
    rcParams['savefig.bbox'] = 'tight'


def plot_target_profiles(data_dir, Z_target=74, Z_ion=1, M_ion=2.0, figsize=(12, 8), plate_filter=None):
    """Main plotting routine."""
    filepath = os.path.join(data_dir, 'target_strike_profile.dat')
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Cannot find {filepath}")

    plates = read_target_data(filepath)
    n_plates = len(plates)
    print(f"Found {n_plates} target plates: {sorted(plates.keys())}")
    for pid, data in plates.items():
        print(f"  Plate {pid}: {data.shape[0]} points, "
              f"strike_dist=[{data[:,1].min():.4f}, {data[:,1].max():.4f}] m")

    # Read strike points if available
    strike_file = os.path.join(data_dir, 'strike_points.dat')
    strike_RZ = []
    if os.path.exists(strike_file):
        strike_RZ = np.loadtxt(strike_file)
        if strike_RZ.ndim == 1:
            strike_RZ = strike_RZ.reshape(1, -1)
        print(f"Read {len(strike_RZ)} strike points from strike_points.dat")
        for i, (sr, sz) in enumerate(strike_RZ):
            print(f"  Strike {i+1}: R={sr:.4f}, Z={sz:.4f}")

    # Column indices: 0:plate_id, 1:strike_dist, 2:R, 3:Z, 4:angle,
    # 5:density, 6:Ti, 7:Te, 8:Vpar, 9:nv, 10:nvT_gam, 11:KparT, 12:psi, 13:sputter

    M_target = 183.84  # W
    colors = ['#0000FF', '#FF0000', '#00FFFF', '#000000']
    labels = ['UOT', 'UIT', 'LIT', 'LOT']

    # Pre-compute shared derived quantities for all plates
    all_sdist_cm = {}
    all_Ti = {}
    all_Te = {}
    all_density = {}
    all_angle_rad = {}
    all_v_sound = {}
    all_Gamma = {}
    all_E_in = {}

    # Determine inner/outer for each plate (compare avg R with x-point R)
    R_xpt = np.mean(strike_RZ[:, 0]) if len(strike_RZ) > 0 else 1.7
    plate_type = {}
    for pid, data in plates.items():
        avg_R = np.mean(data[:, 2])
        plate_type[pid] = 'inner' if avg_R < R_xpt else 'outer'

    # Filter plates if requested
    if plate_filter:
        plates = {pid: data for pid, data in plates.items()
                  if plate_type[pid] == plate_filter}
        if not plates:
            print(f"No plates matching '{args.plate}'")
            return

    for pid, data in plates.items():
        all_sdist_cm[pid] = data[:, 1] * 100  # m -> cm
        all_Ti[pid] = data[:, 6]
        all_Te[pid] = data[:, 7]
        all_density[pid] = np.abs(data[:, 5]) * 10  # 10^19 m^-3
        all_angle_rad[pid] = np.abs(90.0 - data[:, 4]) * np.pi / 180.0
        all_v_sound[pid] = np.sqrt(5.0/3.0 * 2.0 * all_Ti[pid] * 1.602e-19 / (M_ion * 1.67e-27))
        all_Gamma[pid] = all_density[pid] * all_v_sound[pid] * np.sin(all_angle_rad[pid]) * 1e-4  # 10^22 m^-2 s^-1
        all_E_in[pid] = 2.0 * all_Ti[pid] + 3.0 * all_Te[pid]

    def auto_scale(y_dict):
        """Scale y data so axis numbers stay in 0.1-10 range."""
        y_all = np.concatenate(list(y_dict.values()))
        ymax = np.max(np.abs(y_all))
        if ymax <= 0:
            return y_dict, ''
        exp = np.floor(np.log10(ymax))
        if -1 <= exp <= 3:
            return y_dict, ''
        scale = 10.0**(-exp)
        return {k: v * scale for k, v in y_dict.items()}, f' ($\\times 10^{{{int(exp)}}}$)'

    def plot_one(ax, y_data_dict, ylabel, xlim_pad=1.2, ylim_pad=1.5):
        """Plot one quantity for all plates with adaptive axes, MATLAB-style box."""
        y_scaled, scale_suffix = auto_scale(y_data_dict)

        all_x, all_y = [], []
        for i, (pid, data) in enumerate(sorted(plates.items())):
            x = all_sdist_cm[pid]
            y = y_scaled[pid]
            color = colors[i % len(colors)]
            lbl = labels[pid-1] if pid <= len(labels) else f'Plate {pid}'
            ax.plot(x, y, color=color, linewidth=4, label=lbl)
            all_x.append(x)
            all_y.append(y)

        x_all = np.concatenate(all_x)
        y_all = np.concatenate(all_y)
        xmin, xmax = get_xlim(x_all)
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([0, np.max(y_all) * ylim_pad])

        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['top'].set_linewidth(4)
        ax.spines['right'].set_linewidth(4)
        ax.spines['left'].set_linewidth(4)
        ax.spines['bottom'].set_linewidth(4)
        ax.set_box_aspect(1.15/1.45)
        ax.tick_params(top=True, right=True, which='both')

        ax.set_xlabel(r'$R - R_{\mathrm{strike}}$ (cm)')
        ax.set_ylabel(ylabel + scale_suffix)
        ax.legend(loc='best')
        return xmin, xmax

    setup_style()
    # Reserve fixed margins so all figures have identical dimensions
    fig_margin = {'left': 0.073, 'right': 0.927, 'top': 0.95, 'bottom': 0.18}

    # 1. Temperature
    fig1, ax1 = plt.subplots(figsize=figsize)
    plot_one(ax1, all_Ti, r'$T_i$ (eV)', 'Ti')
    fig1.subplots_adjust(**fig_margin)
    fig1.savefig(os.path.join(data_dir, 'target_Ti.png'), dpi=300)

    # 2. Density
    fig2, ax2 = plt.subplots(figsize=figsize)
    plot_one(ax2, all_density, r'$n_e$ ($10^{19}$ m$^{-3}$)', 'density')
    fig2.subplots_adjust(**fig_margin)
    fig2.savefig(os.path.join(data_dir, 'target_density.png'), dpi=300)

    # 3. Particle flux
    fig3, ax3 = plt.subplots(figsize=figsize)
    plot_one(ax3, all_Gamma, r'$\Gamma_\bot$ ($10^{22}$ m$^{-2}$s$^{-1}$)', 'flux')
    fig3.subplots_adjust(**fig_margin)
    fig3.savefig(os.path.join(data_dir, 'target_flux.png'), dpi=300)

    # 4. Heat flux
    H_flux = {pid: all_Gamma[pid] * 1e22 * (8.0 * all_Ti[pid] * 1.602e-19) * 1e-6
              for pid in plates}
    fig4, ax4 = plt.subplots(figsize=figsize)
    plot_one(ax4, H_flux, r'Heat flux (MW/m$^2$)', 'heat flux')
    fig4.subplots_adjust(**fig_margin)
    fig4.savefig(os.path.join(data_dir, 'target_heatflux.png'), dpi=300)

    # 5. Sputtering yield
    yield_phys = {}
    for pid in plates:
        flux = all_Gamma[pid] * 1e22
        yield_phys[pid] = eckstein_sputter_yield(
            Z_ion, Z_target, M_ion, M_target,
            all_E_in[pid], all_angle_rad[pid] * 180/np.pi)
        if Z_target == 6:
            yield_chem = chemical_sputtering(np.maximum(all_Ti[pid], 0.01), all_E_in[pid], flux)
            yield_phys[pid] = yield_phys[pid] + yield_chem

    fig5, ax5 = plt.subplots(figsize=figsize)
    plot_one(ax5, {pid: y * 100 for pid, y in yield_phys.items()},
             r'Sputtering yield ($\times 10^{-2}$)', 'sputtering')
    fig5.subplots_adjust(**fig_margin)
    ax5.yaxis.label.set_size(21)
    fig5.savefig(os.path.join(data_dir, 'target_sputtering.png'), dpi=300)

    # 6. Erosion flux
    erosion = {pid: yield_phys[pid] * all_Gamma[pid] * 1e22 * 1e-21 for pid in plates}
    fig6, ax6 = plt.subplots(figsize=figsize)
    plot_one(ax6, erosion, r'Gross erosion flux ($10^{21}$ m$^{-2}$s$^{-1}$)', 'erosion')
    fig6.subplots_adjust(**fig_margin)
    ax6.yaxis.label.set_size(21)
    fig6.savefig(os.path.join(data_dir, 'target_erosion.png'), dpi=300)

    print("All figures saved.")
    plt.show()


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Plot target plate profiles from jorek2_target2vtkn output.')
    parser.add_argument('data_dir', help='Directory containing target_strike_profile.dat')
    parser.add_argument('--Z_ion', type=int, default=1, help='Ion atomic number (1=H/D/T)')
    parser.add_argument('--Z_target', type=int, default=74, help='Target atomic number (74=W)')
    parser.add_argument('--M_ion', type=float, default=2.0, help='Ion mass in amu (2.0=D)')
    parser.add_argument('--figsize', type=str, default='12,8', help='Figure size W,H in inches')
    parser.add_argument('--plate', type=str, default=None,
                        help='Filter by plate type: inner (R < R_xpoint) or outer (R > R_xpoint)')

    args = parser.parse_args()
    w, h = map(float, args.figsize.split(','))

    plot_target_profiles(args.data_dir, args.Z_target, args.Z_ion, args.M_ion,
                         figsize=(w, h), plate_filter=args.plate)


if __name__ == '__main__':
    main()
