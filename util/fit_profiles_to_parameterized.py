#!/usr/bin/env python3
"""
fit_profiles_to_parameterized.py

Fit JOREK numerical profile data files (jorek_density, jorek_temperature,
jorek_ffprime) to the parameterized analytical form used in model600/601,
and output the corresponding parameter files.

Usage:
    cd /path/to/run/directory
    python /path/to/fit_profiles_to_parameterized.py [--fullmhd DELTA_PSI]

    --fullmhd DELTA_PSI   Output FF' parameters for full-MHD (model750).
                          DELTA_PSI = psi_bnd - psi_axis from the Grad-Shafranov
                          equilibrium (e.g., from eqdsk or fort.10).
                          Without this flag, FF' output is for reduced-MHD.

The script auto-detects jorek_density, jorek_temperature, jorek_temperature_i,
jorek_temperature_e, jorek_ffprime in the current directory and outputs:
    density_parameter     - rho_0, rho_1, rho_coef(1:5)
    temperature_parameter - T_0/Ti_0/Te_0, T_1/Ti_1/Te_1, T_coef/Ti_coef/Te_coef(1:5)
    ffprime_parameter     - FF_0, FF_1, FF_coef(1:9)

The output format is directly pasteable into the JOREK namelist (&in1).

FF' perturbation scaling (see ffprime.f90):
    Reduced-MHD (default):  FF_coef(9)=0, data fitted in normalized psi_n
    Full-MHD (--fullmhd):   FF_coef(9)=1, FF_coef(6) denormalized by delta_psi
                            FF_coef(6)_fullmhd = FF_coef(6)_reduced / delta_psi
"""

import os
import sys
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
import warnings

# Suppress non-critical fitting warnings
warnings.filterwarnings("ignore", category=OptimizeWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ==============================================================================
# Model Functions (mirroring JOREK's density.f90, temperature.f90, ffprime.f90)
# ==============================================================================

def parameterized_density(psi_n, rho_0, rho_1, c1, c2, c3, sig_n, psi_barrier):
    """
    JOREK density parameterization (density.f90).
    prof = (rho_0 - rho_1) * (1 + c1*psi_n + c2*psi_n^2 + c3*psi_n^3)
           * (0.5 - 0.5*tanh((psi_n - psi_barrier)/sig_n))
           + rho_1
    """
    psi_n = np.clip(psi_n, 0.0, 2.0)
    psi_star = np.clip((psi_n - psi_barrier) / max(sig_n, 1e-6), -40.0, 40.0)
    poly = 1.0 + c1 * psi_n + c2 * psi_n**2 + c3 * psi_n**3
    atn = 0.5 - 0.5 * np.tanh(psi_star)
    return (rho_0 - rho_1) * poly * atn + rho_1


def parameterized_temperature(psi_n, T_0, T_1, c1, c2, c3, sig_T, psi_barrier):
    """
    JOREK temperature parameterization (temperature.f90).
    Same functional form as density.
    """
    return parameterized_density(psi_n, T_0, T_1, c1, c2, c3, sig_T, psi_barrier)


def parameterized_ffprime(psi_n, FF_0, FF_1, c1, c2, c3, sig_F, psi_barrier,
                          pert_amp, pert_pos, pert_width):
    """
    JOREK FF' parameterization (ffprime.f90).
    prof = [(FF_0 - FF_1) * (1 + c1*psi_n + c2*psi_n^2 + c3*psi_n^3)
            + pert_amp / cosh((psi_n - pert_pos)/pert_width)^2]
           * (0.5 - 0.5*tanh((psi_n - psi_barrier)/sig_F))
           + FF_1
    """
    psi_n = np.clip(psi_n, 0.0, 2.0)
    psi_star = np.clip((psi_n - psi_barrier) / max(sig_F, 1e-6), -40.0, 40.0)
    poly = 1.0 + c1 * psi_n + c2 * psi_n**2 + c3 * psi_n**3
    pert = pert_amp / np.cosh((psi_n - pert_pos) / max(pert_width, 1e-6))**2
    prof0 = (FF_0 - FF_1) * poly + pert
    atn = 0.5 - 0.5 * np.tanh(psi_star)
    return prof0 * atn + FF_1


# ==============================================================================
# Utility
# ==============================================================================

def auto_detect_delta_psi(work_dir):
    """
    Try to auto-detect delta_psi = psi_bnd - psi_axis in JOREK units.
    Priority:
      1. fort.10 — JOREK equilibrium output (exact JOREK units)
      2. eqdsk.dat — GEQDSK file (SI units, uses simag/sibry as fallback)
    Returns delta_psi or None.
    """
    # --- Try fort.10 (JOREK equilibrium output) ---
    fort10_path = os.path.join(work_dir, 'fort.10')
    if os.path.exists(fort10_path):
        try:
            with open(fort10_path, 'r') as f:
                first_line = f.readline()
            parts = first_line.split()
            # Format: nvar, psi_axis, psi_bnd, psi_xpoint, F0
            # psi_axis and psi_bnd are in JOREK-normalized units
            if len(parts) >= 3:
                psi_axis = float(parts[1])
                psi_bnd  = float(parts[2])
                dpsi = abs(psi_bnd - psi_axis)
                if dpsi > 1e-12:
                    return dpsi
        except Exception:
            pass

    # --- Try eqdsk.dat (GEQDSK, SI units — less precise fallback) ---
    eqdsk_path = os.path.join(work_dir, 'eqdsk.dat')
    if os.path.exists(eqdsk_path):
        try:
            with open(eqdsk_path, 'r') as f:
                # Skip header (line 1), read line 2 (rdim,zdim,rcentr,rleft,zmid)
                f.readline()
                # Line 3: raxis, zaxis, simag, sibry, bcentr
                # Some eqdsk have 2 lines of header; try both patterns
                line3 = f.readline()
            parts = line3.split()
            if len(parts) >= 5:
                simag = float(parts[2])
                sibry = float(parts[3])
                dpsi_si = abs(sibry - simag)
                if dpsi_si > 1e-12:
                    print(f"  Note: using eqdsk SI delta_psi = {dpsi_si:.6g} (Wb).")
                    print(f"  This is approximate; for JOREK units use fort.10 or --delta-psi.")
                    return dpsi_si
        except Exception:
            pass

    return None


def read_profile(filename):
    """Read a two-column ASCII profile file (psi_n, value)."""
    data = np.loadtxt(filename)
    psi_n = data[:, 0]
    values = data[:, 1]
    return psi_n, values


def write_parameter_file(filename, params_dict, description, mode='w'):
    """Write fitted parameters to a file, JOREK namelist-ready."""
    with open(filename, mode) as f:
        if mode == 'w':
            f.write(f"! {description} — fitted from numerical profile\n")
            f.write(f"! Generated by fit_profiles_to_parameterized.py\n\n")
        else:
            f.write(f"\n! --- {description} ---\n")
        for key, value in params_dict.items():
            if isinstance(value, (list, tuple, np.ndarray)):
                vals_str = ", ".join(f"{v:.12g}" for v in value)
                f.write(f"{key} = {vals_str}\n")
            elif isinstance(value, int):
                f.write(f"{key} = {value}\n")
            else:
                f.write(f"{key} = {value:.12g}\n")
    if mode == 'w':
        print(f"  -> {filename}")


def plot_comparison(psi_n, data, fitted, label, filename, psi_n_fine=None, fitted_fine=None):
    """
    Plot original data vs fitted parameterized profile.
    Saves to filename.png.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7),
                                    gridspec_kw={'height_ratios': [3, 1]})

    # --- Top panel: profiles ---
    ax1.plot(psi_n, data, 'b-', linewidth=2.0, label='Original data')
    if psi_n_fine is not None and fitted_fine is not None:
        ax1.plot(psi_n_fine, fitted_fine, 'r--', linewidth=1.8,
                 label='Parameterized fit')
    else:
        ax1.plot(psi_n, fitted, 'r--', linewidth=1.8, label='Parameterized fit')
    ax1.set_ylabel(label)
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    ax1.set_title(f'{label} — Data vs Parameterized Fit')

    # --- Bottom panel: residual ---
    residual = data - fitted
    ax2.plot(psi_n, residual, 'k-', linewidth=1.0)
    ax2.axhline(y=0, color='gray', linestyle=':', linewidth=0.8)
    ax2.set_xlabel(r'$\psi_N$')
    ax2.set_ylabel('Residual')
    ax2.grid(True, alpha=0.3)

    max_res = np.max(np.abs(residual))
    r2 = 1.0 - np.sum(residual**2) / np.sum((data - np.mean(data))**2)
    ax2.set_title(f'Residual  (max |error| = {max_res:.4g},  R² = {r2:.6f})')

    plt.tight_layout()
    plt.savefig(filename + '.png', dpi=150)
    plt.close()
    print(f"  -> {filename}.png")


# ==============================================================================
# Fitting Logic
# ==============================================================================

def fit_density(psi_n, data):
    """
    Fit density profile.
    7 parameters: rho_0, rho_1, c1, c2, c3, sig_n, psi_barrier
    """
    rho_axis = data[0]   # value at psi_n ≈ 0
    rho_sol  = data[-1]  # value at psi_n ≈ 1.0

    # Robust initial guess
    p0 = [rho_axis, rho_sol, -0.5, 0.3, 0.0, 0.05, 1.0]
    bounds = (
        [rho_axis * 0.5,  0.0,         -5.0, -3.0, -3.0, 0.005, 0.85],
        [rho_axis * 1.5,  rho_sol*2.0,  5.0,  3.0,  3.0, 0.3,   1.2]
    )

    popt, pcov = curve_fit(parameterized_density, psi_n, data, p0=p0,
                           bounds=bounds, maxfev=50000, method='trf')
    rho_0, rho_1, c1, c2, c3, sig_n, psi_barrier = popt

    # JOREK's rho_coef array: (1:5) = c1, c2, c3, sig_n, psi_barrier
    # (indices 6:10 are unused)
    params = {
        'rho_0':     rho_0,
        'rho_1':     rho_1,
        'rho_coef(1)': c1,
        'rho_coef(2)': c2,
        'rho_coef(3)': c3,
        'rho_coef(4)': sig_n,
        'rho_coef(5)': psi_barrier,
    }

    # Compute fit quality
    fitted = parameterized_density(psi_n, *popt)
    r2 = 1.0 - np.sum((data - fitted)**2) / np.sum((data - np.mean(data))**2)
    print(f"  rho_0={rho_0:.6g}, rho_1={rho_1:.6g}, "
          f"c=[{c1:.4f},{c2:.4f},{c3:.4f}], "
          f"sig_n={sig_n:.4f}, barrier={psi_barrier:.4f}  R²={r2:.6f}")

    return params


def fit_temperature(psi_n, data, prefix='T'):
    """
    Fit temperature profile (single T, Ti, or Te).
    7 parameters: T_0, T_1, c1, c2, c3, sig_T, psi_barrier
    prefix = 'T' | 'Ti' | 'Te' controls output parameter names.
    """
    T_axis = data[0]
    T_sol  = data[-1]

    p0 = [T_axis, T_sol, -0.5, 0.3, 0.0, 0.05, 1.0]
    bounds = (
        [T_axis * 0.5,  0.0,        -5.0, -3.0, -3.0, 0.005, 0.85],
        [T_axis * 1.5,  T_sol*2.0,   5.0,  3.0,  3.0, 0.3,   1.2]
    )

    popt, pcov = curve_fit(parameterized_temperature, psi_n, data, p0=p0,
                           bounds=bounds, maxfev=50000, method='trf')
    T_0, T_1, c1, c2, c3, sig_T, psi_barrier = popt

    params = {
        f'{prefix}_0':       T_0,
        f'{prefix}_1':       T_1,
        f'{prefix}_coef(1)': c1,
        f'{prefix}_coef(2)': c2,
        f'{prefix}_coef(3)': c3,
        f'{prefix}_coef(4)': sig_T,
        f'{prefix}_coef(5)': psi_barrier,
    }

    fitted = parameterized_temperature(psi_n, *popt)
    r2 = 1.0 - np.sum((data - fitted)**2) / np.sum((data - np.mean(data))**2)
    print(f"  {prefix}_0={T_0:.6g}, {prefix}_1={T_1:.6g}, "
          f"c=[{c1:.4f},{c2:.4f},{c3:.4f}], "
          f"sig_T={sig_T:.4f}, barrier={psi_barrier:.4f}  R²={r2:.6f}")

    return params


def fit_ffprime(psi_n, data, delta_psi=None, fullmhd=False):
    """
    Fit FF' profile.
    10 parameters: FF_0, FF_1, c1, c2, c3, sig_F, psi_barrier,
                   pert_amp, pert_pos, pert_width

    delta_psi = None   → reduced-MHD w/o delta_psi (warning, may need correction)
    delta_psi = value, fullmhd=False → reduced-MHD with corrected FF_coef(6)
    delta_psi = value, fullmhd=True  → full-MHD, FF_coef(9)=1

    JOREK formula (ffprime.f90):
      d_pert = FF_coef(6)/cosh(...)^2 / (2*FF_coef(8)) / delta_psi * no_delta_psi
      no_delta_psi = 1 (FF_coef(9)=0) or delta_psi (FF_coef(9)=1)

    Conversion:
      reduced-MHD: FF_coef(6) = pert_amp * 2 * pert_width * delta_psi
      full-MHD:    FF_coef(6) = pert_amp * 2 * pert_width
    """
    FF_axis = data[0]
    FF_sol  = data[-1]

    # Estimate perturbation: residual after smooth trend in edge region
    edge_mask = psi_n > 0.8
    if np.any(edge_mask):
        pert_guess = np.max(np.abs(data[edge_mask] - np.mean(data[edge_mask])))
    else:
        pert_guess = 0.001

    p0 = [FF_axis, FF_sol, -0.5, 0.3, 0.0, 0.05, 1.0,
          pert_guess, 0.98, 0.05]
    bounds = (
        [FF_axis * 2.0,  FF_sol * 2.0, -5.0, -3.0, -3.0, 0.005, 0.85,
         0.0,            0.88,          0.005],
        [FF_axis * 0.5,  FF_sol * 0.5,  5.0,  3.0,  3.0, 0.3,   1.2,
         abs(FF_axis)*5, 1.05,          0.3]
    )
    # Note: FF' is negative, so FF_0 < 0 and FF_1 < 0 typically.
    # Adjust bound ordering accordingly
    FF_0_lower = min(FF_axis * 2.0, FF_axis * 0.5)
    FF_0_upper = max(FF_axis * 2.0, FF_axis * 0.5)
    FF_1_lower = min(FF_sol * 2.0, FF_sol * 0.5)
    FF_1_upper = max(FF_sol * 2.0, FF_sol * 0.5)

    bounds = (
        [FF_0_lower,  FF_1_lower, -5.0, -3.0, -3.0, 0.005, 0.85,
         0.0,          0.88,       0.005],
        [FF_0_upper,  FF_1_upper,  5.0,  3.0,  3.0, 0.3,   1.2,
         abs(FF_axis)*5, 1.05,     0.3]
    )

    try:
        popt, pcov = curve_fit(parameterized_ffprime, psi_n, data, p0=p0,
                               bounds=bounds, maxfev=100000, method='trf')
    except Exception as e:
        print(f"  WARNING: Full FF' fit failed ({e}), trying without perturbation...")
        # Fallback: fit without perturbation (set pert_amp=0)
        def ffprime_no_pert(psi_n, FF_0, FF_1, c1, c2, c3, sig_F, psi_barrier):
            return parameterized_ffprime(psi_n, FF_0, FF_1, c1, c2, c3,
                                         sig_F, psi_barrier, 0.0, 1.0, 0.1)
        p0_simple = [FF_axis, FF_sol, -0.5, 0.3, 0.0, 0.05, 1.0]
        bounds_simple = (
            [FF_0_lower,  FF_1_lower, -5.0, -3.0, -3.0, 0.005, 0.85],
            [FF_0_upper,  FF_1_upper,  5.0,  3.0,  3.0, 0.3,   1.2]
        )
        popt, pcov = curve_fit(ffprime_no_pert, psi_n, data, p0=p0_simple,
                               bounds=bounds_simple, maxfev=50000, method='trf')
        FF_0, FF_1, c1, c2, c3, sig_F, psi_barrier = popt
        pert_amp, pert_pos, pert_width = 0.0, 1.0, 0.1
    else:
        FF_0, FF_1, c1, c2, c3, sig_F, psi_barrier, pert_amp, pert_pos, pert_width = popt

    # --- Convert to JOREK's FF_coef convention ---
    # JOREK formula (ffprime.f90 line 77):
    #   d_pert = FF_coef(6) / cosh(...)^2 / (2*FF_coef(8)) / delta_psi * no_delta_psi
    #   where no_delta_psi = 1           if FF_coef(9)=0 (reduced-MHD)
    #         no_delta_psi = delta_psi   if FF_coef(9)=1 (full-MHD)
    #
    # We fitted: pert = pert_amp / cosh((psi_n - pert_pos)/pert_width)^2
    # So:  pert_amp = FF_coef(6) / (2*pert_width) / delta_psi * no_delta_psi
    #  =>  FF_coef(6) = pert_amp * 2 * pert_width * delta_psi / no_delta_psi
    #
    # Reduced-MHD (FF_coef(9)=0, no_delta_psi=1):
    #   FF_coef(6)_reduced = pert_amp * 2 * pert_width * delta_psi   ← needs delta_psi
    # Full-MHD  (FF_coef(9)=1, no_delta_psi=delta_psi):
    #   FF_coef(6)_fullmhd = pert_amp * 2 * pert_width               ← delta_psi cancels
    # Conversion: FF_coef(6)_fullmhd = FF_coef(6)_reduced / delta_psi
    pert_amp_jorek = pert_amp * 2.0 * pert_width

    if fullmhd:
        if delta_psi is None or delta_psi <= 0:
            print("  ERROR: --fullmhd requires a positive delta_psi value.")
            sys.exit(1)
        ff_coef_6 = pert_amp_jorek
        ff_coef_9 = 1
        mode_note = f"full-MHD (FF_coef(9)=1, denormalized by delta_psi={delta_psi:.6g})"
    elif delta_psi is not None:
        if delta_psi <= 0:
            print("  ERROR: delta_psi must be positive.")
            sys.exit(1)
        ff_coef_6 = pert_amp_jorek * delta_psi
        ff_coef_9 = 0
        mode_note = f"reduced-MHD (FF_coef(9)=0, delta_psi={delta_psi:.6g})"
    else:
        ff_coef_6 = pert_amp_jorek
        ff_coef_9 = 0
        mode_note = "reduced-MHD (WARNING: delta_psi unknown, FF_coef(6) may need correction)"
        print(f"  !! WARNING: delta_psi not provided.")
        print(f"  !! FF_coef(6) = {ff_coef_6:.6g} assumes delta_psi=1.0.")
        print(f"  !! For correct reduced-MHD results, multiply FF_coef(6) by actual")
        print(f"  !! delta_psi (= psi_bnd - psi_axis). Or use --delta-psi VALUE.")
        print(f"  !! For full-MHD (model750), use --fullmhd DELTA_PSI instead.")

    params = {
        'FF_0':       FF_0,
        'FF_1':       FF_1,
        'FF_coef(1)': c1,
        'FF_coef(2)': c2,
        'FF_coef(3)': c3,
        'FF_coef(4)': sig_F,
        'FF_coef(5)': psi_barrier,
        'FF_coef(6)': ff_coef_6,
        'FF_coef(7)': pert_pos,
        'FF_coef(8)': pert_width,
        'FF_coef(9)': ff_coef_9,   # 0=normalized, 1=denormalized
    }

    fitted = parameterized_ffprime(psi_n, *popt)
    r2 = 1.0 - np.sum((data - fitted)**2) / np.sum((data - np.mean(data))**2)
    print(f"  FF_0={FF_0:.6g}, FF_1={FF_1:.6g}, "
          f"c=[{c1:.4f},{c2:.4f},{c3:.4f}], "
          f"sig_F={sig_F:.4f}, barrier={psi_barrier:.4f}, "
          f"pert=[{pert_amp:.4g},{pert_pos:.4f},{pert_width:.4f}]  R²={r2:.6f}")
    print(f"  => {mode_note}")

    return params


# ==============================================================================
# Main
# ==============================================================================

def main():
    # --- Parse command-line arguments ---
    import argparse
    parser = argparse.ArgumentParser(
        description="Fit JOREK numerical profiles to parameterized analytical form.",
        epilog="Output files are directly pasteable into JOREK namelist (&in1)."
    )
    parser.add_argument('--delta-psi', type=float, default=None,
                        help="psi_bnd - psi_axis from Grad-Shafranov equilibrium "
                             "(JOREK units). Auto-detected if not provided.")
    parser.add_argument('--fullmhd', action='store_true', default=None,
                        help="Output FF' for full-MHD (model750): FF_coef(9)=1, "
                             "FF_coef(6) denormalized. Uses auto-detected delta_psi.")
    parser.add_argument('--fullmhd-dpsi', type=float, metavar='DELTA_PSI', default=None,
                        help="Same as --fullmhd but specify delta_psi manually.")
    args = parser.parse_args()

    if args.fullmhd_dpsi is not None:
        delta_psi = args.fullmhd_dpsi
        fullmhd_mode = True
    elif args.fullmhd:
        delta_psi = None  # auto-detect
        fullmhd_mode = True
    elif args.delta_psi is not None:
        delta_psi = args.delta_psi
        fullmhd_mode = False
    else:
        delta_psi = None  # auto-detect
        fullmhd_mode = False

    # Determine working directory
    work_dir = os.getcwd()
    print(f"Working directory: {work_dir}")

    # --- Auto-detect delta_psi if not provided ---
    if delta_psi is None:
        detected = auto_detect_delta_psi(work_dir)
        if detected is not None:
            delta_psi = detected
            print(f"Auto-detected delta_psi = {delta_psi:.6g} (JOREK units)")
        else:
            print("Could not auto-detect delta_psi from fort.10 or eqdsk.dat.")

    if fullmhd_mode:
        if delta_psi is not None:
            print(f"Full-MHD mode: FF_coef(9)=1, delta_psi={delta_psi:.6g}")
        else:
            print("Full-MHD mode: FF_coef(9)=1, delta_psi UNKNOWN (FF_coef(6) may be wrong!)")
    elif delta_psi is not None:
        print(f"Reduced-MHD mode with delta_psi={delta_psi:.6g}")
    else:
        print("Reduced-MHD mode (delta_psi unknown — may need manual FF_coef(6) correction)")
    print(f"Looking for profile data files...\n")

    # File definitions
    files_config = {
        'jorek_density':     ('density_parameter',     fit_density,     'Density'),
        'jorek_ffprime':     ('ffprime_parameter',     fit_ffprime,     'FFprime'),
    }

    # Temperature files (all write to a single temperature_parameter file)
    temp_files = [
        ('jorek_temperature',   fit_temperature, 'T',  'Temperature'),
        ('jorek_temperature_i', fit_temperature, 'Ti', 'Temperature (ion)'),
        ('jorek_temperature_e', fit_temperature, 'Te', 'Temperature (electron)'),
    ]

    results = {}
    temp_params_list = []  # collect (desc, params) for each temp file found

    for data_file, (out_file, fit_func, desc) in files_config.items():
        data_path = os.path.join(work_dir, data_file)
        if not os.path.exists(data_path):
            print(f"  [{desc}] {data_file} NOT FOUND — skipping.")
            continue

        print(f"  [{desc}] Fitting {data_file}...")
        try:
            psi_n, data = read_profile(data_path)
            if data_file == 'jorek_ffprime':
                params = fit_func(psi_n, data, delta_psi=delta_psi,
                                  fullmhd=fullmhd_mode)
                # Compute fitted values for plotting
                fitted = parameterized_ffprime(psi_n, params['FF_0'], params['FF_1'],
                    params['FF_coef(1)'], params['FF_coef(2)'], params['FF_coef(3)'],
                    params['FF_coef(4)'], params['FF_coef(5)'],
                    params['FF_coef(6)'], params['FF_coef(7)'], params['FF_coef(8)'])
                # Use FF_coef(9) to determine the effective pert_amp for plotting
                # (the fitted function always uses normalized form, so no correction needed for plot)
                plot_comparison(psi_n, data, fitted, 'FF\'', out_file)
            elif data_file == 'jorek_density':
                params = fit_func(psi_n, data)
                fitted = parameterized_density(psi_n, params['rho_0'], params['rho_1'],
                    params['rho_coef(1)'], params['rho_coef(2)'], params['rho_coef(3)'],
                    params['rho_coef(4)'], params['rho_coef(5)'])
                plot_comparison(psi_n, data, fitted, 'Density', out_file)
            else:
                params = fit_func(psi_n, data)
            write_parameter_file(os.path.join(work_dir, out_file), params, desc)
            results[desc] = params
        except Exception as e:
            print(f"  [{desc}] ERROR: {e}")

    # --- Temperature files: all go into temperature_parameter ---
    first_temp = True
    for data_file, fit_func, prefix, desc in temp_files:
        data_path = os.path.join(work_dir, data_file)
        if not os.path.exists(data_path):
            continue

        print(f"  [{desc}] Fitting {data_file}...")
        try:
            psi_n, data = read_profile(data_path)
            params = fit_func(psi_n, data, prefix=prefix)
            param_0 = params[f'{prefix}_0']
            param_1 = params[f'{prefix}_1']
            c1, c2, c3 = params[f'{prefix}_coef(1)'], params[f'{prefix}_coef(2)'], params[f'{prefix}_coef(3)']
            sig, barrier = params[f'{prefix}_coef(4)'], params[f'{prefix}_coef(5)']
            fitted = parameterized_temperature(psi_n, param_0, param_1, c1, c2, c3, sig, barrier)
            plot_comparison(psi_n, data, fitted, f'{prefix} (Temperature)',
                            f'temperature_{prefix.lower()}_parameter')
            mode = 'w' if first_temp else 'a'
            write_parameter_file(os.path.join(work_dir, 'temperature_parameter'),
                                params, desc, mode=mode)
            first_temp = False
            results[desc] = params
        except Exception as e:
            print(f"  [{desc}] ERROR: {e}")

    # Summary
    print(f"\n{'='*60}")
    if results:
        print("  All done. Output files written:")
        if 'Density' in results:
            print("    density_parameter")
        if any('Temperature' in k for k in results):
            print("    temperature_parameter")
        if 'FFprime' in results:
            print("    ffprime_parameter")
        print(f"\n  To use these in JOREK, set in your namelist (&in1):")
        print(f"    rho_file = 'none'")
        if any('Temperature (ion)' in k for k in results) or any('Temperature (electron)' in k for k in results):
            print(f"    Ti_file  = 'none'")
            print(f"    Te_file  = 'none'")
        else:
            print(f"    T_file   = 'none'")
        print(f"    ffprime_file = 'none'")
        if fullmhd_mode:
            print(f"\n  Full-MHD mode: FF_coef(9)=1, FF_coef(6) already denormalized.")
            print(f"  Copy the ffprime_parameter values directly to your namelist.")
        elif delta_psi is None:
            print(f"\n  Reduced-MHD: FF_coef(9)=0. If FF' profile looks wrong,")
            print(f"  rerun with --delta-psi VALUE to correct FF_coef(6).")
        print(f"  And copy the parameter values from the output files.")
    else:
        print("  No profiles were fitted. Check that the data files exist.")
        sys.exit(1)


if __name__ == '__main__':
    main()
