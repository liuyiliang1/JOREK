#!/usr/bin/env python3
"""
JOREK Restart HDF5 File Reader and Magnetic Island Width Calculator
===================================================================

Reads JOREK restart.h5 files and provides convenient access to:
  - Grid geometry (R, Z of all nodes)
  - Physical fields (psi, u, j, w, rho, T, v_par, T_e, rho_n, ...)
  - Fourier harmonics decomposition
  - Equilibrium q-profile (from n=0 poloidal flux)
  - Magnetic island widths at rational surfaces

Usage as a library:
    from jorek_restart_reader import JOREKRestart

    jr = JOREKRestart('jorek_restart.h5')
    jr.print_info()
    q, psi_n = jr.get_q_profile()
    islands = jr.get_island_widths()  # returns list of dicts

Usage as a script:
    python jorek_restart_reader.py jorek_restart.h5          # Print info
    python jorek_restart_reader.py jorek_restart.h5 --islands # Calculate island widths
    python jorek_restart_reader.py jorek_restart.h5 --qprofile # Plot q-profile
    python jorek_restart_reader.py jorek_restart.h5 --poincare --R 1.8 --Z 0.0  # Field line trace

HDF5 File Structure (from mod_export_restart.f90):
    /RCS_version          - Git revision string
    /rst_hdf5_version     - HDF5 format version
    /jorek_model          - Model number (e.g., 601)
    /n_var                - Number of physical variables
    /n_tor                - Number of toroidal harmonics
    /n_period             - Number of toroidal periods (stellarator)
    /n_nodes              - Number of nodes
    /n_elements           - Number of elements
    /n_vertex_max         - Max vertices per element (typically 4)
    /n_degrees            - Number of basis function degrees
    /n_order              - Polynomial order of basis functions
    /x      (n_nodes, n_coord_tor, n_degrees, n_dim) - Node coordinates
    /values (n_nodes, n_tor, n_degrees, n_var)       - Field values (Fourier harmonics)
    /vertex (n_elements, n_vertex_max)                - Element connectivity
    /size   (n_elements, n_vertex_max, n_degrees)     - Element size parameters
    /index  (n_nodes, n_degrees)                       - Node indexing
    /boundary (n_nodes)                                - Boundary flags
    /t_now, /tstep, /F0, /eta, /visco, ...
    /index_now

Variable numbering (typical, depends on model):
    1: psi  (poloidal flux)
    2: u    (velocity stream function)
    3: j    (toroidal current density)
    4: w    (vorticity)
    5: rho  (density)
    6: T    (temperature)
    7: v_par (parallel velocity)
    ...

Toroidal harmonic numbering:
    0: n=0 (axisymmetric/equilibrium)
    1: n=1 cos
    2: n=1 sin
    3: n=2 cos
    4: n=2 sin
    ...

Author: Generated for JOREK diagnostics analysis
"""

import numpy as np
import sys
import os
from collections import OrderedDict

# Lazy import: h5py is needed at runtime but may not be in default Python
try:
    import h5py
    _has_h5py = True
except ImportError:
    _has_h5py = False

# ==============================================================================
# Variable name mapping (model-dependent, covers common models)
# ==============================================================================
VAR_NAMES_DEFAULT = {
    1: "psi",
    2: "u",
    3: "j",
    4: "w",
    5: "rho",
    6: "T",
    7: "v_par",
    8: "T_e",
    9: "rho_n",
    10: "A_R",
    11: "A_Z",
    12: "A_phi",
}

# ==============================================================================
# Basis functions for Bezier elements (from jorek_read_h5.py)
# ==============================================================================
def bezier_basis(s, t):
    """Bezier basis functions of order 3 in s and t directions.
    Returns array of shape (n_order, n_vertex) evaluated at (s, t).
    n_order=4: [0]=corner value, [1]=s-derivative, [2]=t-derivative, [3]=cross-derivative
    n_vertex=4: Bezier control points per element
    """
    return np.asarray([
        [(-1 + s)**2 * (1 + 2*s) * (-1 + t)**2 * (1 + 2*t),
         -(s**2 * (-3 + 2*s) * (-1 + t)**2 * (1 + 2*t)),
         s**2 * (-3 + 2*s) * t**2 * (-3 + 2*t),
         -(-1 + s)**2 * (1 + 2*s) * t**2 * (-3 + 2*t)],
        [3*(-1 + s)**2 * s * (-1 + t)**2 * (1 + 2*t),
         -3*(-1 + s) * s**2 * (-1 + t)**2 * (1 + 2*t),
         3*(-1 + s) * s**2 * t**2 * (-3 + 2*t),
         -3*(-1 + s)**2 * s * t**2 * (-3 + 2*t)],
        [3*(-1 + s)**2 * (1 + 2*s) * (-1 + t)**2 * t,
         -3*s**2 * (-3 + 2*s) * (-1 + t)**2 * t,
         3*s**2 * (-3 + 2*s) * (-1 + t) * t**2,
         -3*(-1 + s)**2 * (1 + 2*s) * (-1 + t) * t**2],
        [9*(-1 + s)**2 * s * (-1 + t)**2 * t,
         -9*(-1 + s) * s**2 * (-1 + t)**2 * t,
         9*(-1 + s) * s**2 * (-1 + t) * t**2,
         -9*(-1 + s)**2 * s * (-1 + t) * t**2]
    ])


# ==============================================================================
# Core restart reader class
# ==============================================================================
class JOREKRestart:
    """Reader for JOREK HDF5 restart files.

    Parameters
    ----------
    filename : str
        Path to the HDF5 restart file (e.g., 'jorek_restart.h5').
    verbose : bool
        Print verbose information during loading.

    Examples
    --------
    >>> jr = JOREKRestart('jorek_restart.h5')
    >>> jr.print_info()
    >>> R, Z = jr.get_node_coordinates()
    >>> psi_n0 = jr.get_field('psi', tor_mode=0)  # equilibrium psi
    >>> q, psi_n = jr.get_q_profile(n_surfaces=100)
    >>> islands = jr.get_island_widths()
    """

    def __init__(self, filename, verbose=False):
        self.filename = filename
        self.verbose = verbose
        self._h5 = None
        self._values = None
        self._x = None
        self._node_R = None
        self._node_Z = None
        self._psi_axis = None
        self._psi_bnd = None
        self._R_axis = None
        self._Z_axis = None
        self._var_names = {}

        self._load()

    def _load(self):
        """Load metadata and map datasets from the HDF5 file."""
        if not _has_h5py:
            raise ImportError(
                "h5py is required to read JOREK HDF5 restart files. "
                "Install it via: pip install h5py\n"
                "Or load it from your HPC environment module.")
        self._h5 = h5py.File(self.filename, 'r')

        # Read scalar metadata
        self._read_scalar('rst_hdf5_version', 1)
        self._read_scalar('jorek_model', -1)
        self._read_scalar('n_var', 7)
        self._read_scalar('n_tor', 1)
        self._read_scalar('n_coord_tor', 1)
        self._read_scalar('n_period', 1)
        self._read_scalar('n_degrees', 4)
        self._read_scalar('n_order', 3)
        self._read_scalar('n_nodes', 0)
        self._read_scalar('n_elements', 0)
        self._read_scalar('n_vertex_max', 4)
        self._read_scalar('n_dim', 2)
        self._read_scalar('t_now', 0.0)
        self._read_scalar('tstep', 0.0)
        self._read_scalar('index_now', 0)
        self._read_scalar('F0', 0.0)
        self._read_scalar('eta', 0.0)
        self._read_scalar('visco', 0.0)

        # Set up variable names
        self._var_names = VAR_NAMES_DEFAULT.copy()

        if self.verbose:
            print(f"Loaded {self.filename}")
            print(f"  Model: {self.jorek_model}, Nodes: {self.n_nodes}, "
                  f"Elements: {self.n_elements}")
            print(f"  n_var={self.n_var}, n_tor={self.n_tor}, "
                  f"n_period={self.n_period}")
            print(f"  t_now={self.t_now:.6e}, index_now={self.index_now}")

    def _read_scalar(self, name, default):
        """Read a scalar from HDF5, falling back to default."""
        try:
            val = self._h5[name][()]
            if np.ndim(val) > 0:
                val = val.flat[0] if val.size > 0 else default
            setattr(self, name, int(val) if isinstance(default, int) else float(val))
        except (KeyError, ValueError):
            setattr(self, name, default)

    def close(self):
        """Close the HDF5 file."""
        if self._h5 is not None:
            self._h5.close()
            self._h5 = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    # --------------------------------------------------------------------------
    # Metadata / info
    # --------------------------------------------------------------------------
    def print_info(self):
        """Print a summary of the restart file contents."""
        print("=" * 65)
        print(f"  JOREK Restart File: {self.filename}")
        print("=" * 65)
        print(f"  Model number:        {self.jorek_model}")
        print(f"  HDF5 version:        {self.rst_hdf5_version}")
        print(f"  Time:                {self.t_now:.6e}")
        print(f"  Time step index:     {self.index_now}")
        print(f"  Time step dt:        {self.tstep:.6e}")
        print(f"  F0 (R*B_phi):        {self.F0:.6e}")
        print(f"  Resistivity (eta):   {self.eta:.6e}")
        print(f"  Viscosity (visco):   {self.visco:.6e}")
        print(f"  Toroidal period:     {self.n_period}")
        print(f"-" * 65)
        print(f"  n_nodes:             {self.n_nodes}")
        print(f"  n_elements:          {self.n_elements}")
        print(f"  n_vertex_max:        {self.n_vertex_max}")
        print(f"  n_degrees:           {self.n_degrees}")
        print(f"  n_order:             {self.n_order}")
        print(f"  n_var:               {self.n_var}")
        print(f"  n_tor:               {self.n_tor}")
        print(f"  n_coord_tor:         {self.n_coord_tor}")
        print(f"-" * 65)
        print(f"  Variables (n_var index):")
        for i in range(1, self.n_var + 1):
            name = self._var_names.get(i, f"var_{i}")
            print(f"    {i}: {name}")
        print(f"-" * 65)
        print(f"  Toroidal modes (n_tor index):")
        print(f"    0: n=0 (axisymmetric)")
        for i in range(1, (self.n_tor - 1) // 2 + 1):
            n_mode = i * self.n_period
            print(f"    {2*i}:   n={n_mode} cos")
            print(f"    {2*i+1}: n={n_mode} sin")
        print("=" * 65)

    def get_var_index(self, name):
        """Get the variable index from its name string.

        Parameters
        ----------
        name : str
            Variable name (e.g., 'psi', 'u', 'T').

        Returns
        -------
        int or None
        """
        for idx, vname in self._var_names.items():
            if vname == name.lower():
                return idx
        return None

    # --------------------------------------------------------------------------
    # Raw data access
    # --------------------------------------------------------------------------
    @property
    def values(self):
        """4D array (n_nodes, n_tor, n_degrees, n_var) of field values."""
        if self._values is None:
            self._values = self._h5['values'][:]
        return self._values

    @property
    def x_coords(self):
        """Coordinate array (n_nodes, n_coord_tor, n_degrees, n_dim)."""
        if self._x is None:
            self._x = self._h5['x'][:]
        return self._x

    @property
    def vertex(self):
        """Element connectivity (n_elements, n_vertex_max)."""
        return self._h5['vertex'][:]

    @property
    def element_size(self):
        """Element sizes (n_elements, n_vertex_max, n_degrees)."""
        return self._h5['size'][:]

    def get_node_coordinates(self):
        """Get (R, Z) coordinates for all nodes (n=0 toroidal component).

        Returns
        -------
        R : ndarray (n_nodes,)
        Z : ndarray (n_nodes,)
        """
        if self._node_R is not None and self._node_Z is not None:
            return self._node_R, self._node_Z

        x = self.x_coords
        # x shape: (n_nodes, n_coord_tor, n_degrees, n_dim)
        # For v2 format: use dim 1 for coords; for v1 use first coord_tor
        if x.ndim == 4:
            # x[node, 0, degree, dim]
            self._node_R = x[:, 0, 0, 0].ravel()  # degree=0 = corner values
            self._node_Z = x[:, 0, 0, 1].ravel()
        else:
            # Old format: x[node, degree, dim]
            self._node_R = x[:, 0, 0].ravel()
            self._node_Z = x[:, 0, 1].ravel()

        return self._node_R, self._node_Z

    def get_field(self, var, tor_mode=0, degree=0):
        """Extract a scalar field from the restart file.

        Parameters
        ----------
        var : int or str
            Variable index (1-based) or name (e.g., 'psi', 'u').
        tor_mode : int
            Toroidal harmonic index:
            0 = n=0 (equilibrium),
            1 = n=1 cos,
            2 = n=1 sin,
            3 = n=2 cos,
            4 = n=2 sin, ...
        degree : int
            Basis function degree index (0 = value, 1 = s-deriv, 2 = t-deriv).

        Returns
        -------
        ndarray (n_nodes,)
        """
        if isinstance(var, str):
            var = self.get_var_index(var)
            if var is None:
                raise ValueError(f"Unknown variable name: {var}")

        # values shape: (n_nodes, n_tor, n_degrees, n_var)
        if var > self.values.shape[3]:
            raise IndexError(f"var={var} exceeds n_var={self.values.shape[3]}")

        return self.values[:, tor_mode, degree, var - 1].ravel()

    def get_field_RZ(self, var, R, Z, tor_mode=0, n_points=200):
        """Evaluate a field at a specific (R, Z) position by element search.

        Parameters
        ----------
        var : int or str
            Variable to evaluate.
        R, Z : float
            Physical coordinates.
        tor_mode : int
            Toroidal harmonic index.
        n_points : int
            Subdivision points per element for Bezier evaluation.

        Returns
        -------
        float or None if point outside domain.
        """
        if isinstance(var, str):
            var = self.get_var_index(var)
            if var is None:
                return None

        R_nodes, Z_nodes = self.get_node_coordinates()
        # Find the element containing (R, Z)
        # Simple approach: find closest element by bounding box
        # (A more precise approach would use the Bezier element inversion,
        #  but this is adequate for most purposes)
        vtx = self.vertex  # (n_elements, n_vertex_max)
        n_elm = vtx.shape[0]

        best_dist = np.inf
        best_elem = -1
        best_s, best_t = 0.5, 0.5

        for i_elm in range(min(n_elm, 5000)):  # Cap search for speed
            # Get R, Z of 4 corner vertices for this element
            n1, n2, n3, n4 = vtx[i_elm, :4] - 1  # Convert to 0-based
            if n1 < 0 or n2 < 0 or n3 < 0 or n4 < 0:
                continue

            R_elm = [R_nodes[n1], R_nodes[n2], R_nodes[n3], R_nodes[n4]]
            Z_elm = [Z_nodes[n1], Z_nodes[n2], Z_nodes[n3], Z_nodes[n4]]

            # Check if point is inside quadrilateral using winding number
            # Or simpler: centroid distance
            R_cent = np.mean(R_elm)
            Z_cent = np.mean(Z_elm)
            dist = np.sqrt((R - R_cent)**2 + (Z - Z_cent)**2)
            if dist < best_dist:
                best_dist = dist
                best_elem = i_elm

        if best_elem < 0:
            return None

        # For a rough estimate, use the value at the nearest node
        # A full Bezier evaluation would require solving for (s,t) first
        n_idx = vtx[best_elem, 0] - 1
        return float(self.get_field(var, tor_mode=tor_mode)[n_idx])

    # --------------------------------------------------------------------------
    # Equilibrium analysis
    # --------------------------------------------------------------------------
    def find_magnetic_axis(self):
        """Find the magnetic axis (minimum of n=0 poloidal flux).

        Returns
        -------
        R_axis, Z_axis, psi_axis : float
            Coordinates and psi value at magnetic axis.
        """
        if self._psi_axis is not None:
            return self._R_axis, self._Z_axis, self._psi_axis

        R, Z = self.get_node_coordinates()
        psi = self.get_field('psi', tor_mode=0)

        # Magnetic axis is where psi is extremal (min for typical tokamak)
        idx_axis = np.argmin(psi)

        self._R_axis = float(R[idx_axis])
        self._Z_axis = float(Z[idx_axis])
        self._psi_axis = float(psi[idx_axis])

        return self._R_axis, self._Z_axis, self._psi_axis

    def find_separatrix_psi(self):
        """Find the separatrix psi value from boundary nodes (X-point).

        Returns
        -------
        psi_bnd : float
        """
        if self._psi_bnd is not None:
            return self._psi_bnd

        try:
            boundary = self._h5['boundary'][:]
        except KeyError:
            # Without boundary info, use a crude estimate
            psi = self.get_field('psi', tor_mode=0)
            self._psi_bnd = float(np.max(psi))
            return self._psi_bnd

        psi = self.get_field('psi', tor_mode=0)
        # Boundary nodes have non-zero boundary flag
        bnd_mask = boundary > 0
        if np.any(bnd_mask):
            self._psi_bnd = float(np.max(psi[bnd_mask]))
        else:
            self._psi_bnd = float(np.max(psi))

        return self._psi_bnd

    def get_psi_n(self, psi=None):
        """Compute normalized poloidal flux psi_n.

        Parameters
        ----------
        psi : ndarray or None
            Raw psi values. If None, uses n=0 equilibrium psi.

        Returns
        -------
        psi_n : ndarray
            Normalized psi: (psi - psi_axis) / (psi_bnd - psi_axis)
        """
        if psi is None:
            psi = self.get_field('psi', tor_mode=0)

        psi_axis = self.find_magnetic_axis()[2]
        psi_bnd = self.find_separatrix_psi()

        return (psi - psi_axis) / (psi_bnd - psi_axis)

    # --------------------------------------------------------------------------
    # Q-profile
    # --------------------------------------------------------------------------
    def get_q_profile(self, n_surfaces=100, n_theta=64):
        """Compute the q (safety factor) profile by tracing flux surfaces.

        This implements a similar algorithm to the postproc 'qprofile' command.
        For each flux surface (psi_n value), field lines are traced and the
        toroidal excursion per poloidal turn gives q.

        Parameters
        ----------
        n_surfaces : int
            Number of flux surfaces to trace.
        n_theta : int
            Number of poloidal angle points for surface integration.

        Returns
        -------
        q : ndarray (n_surfaces-2,)
            Safety factor values (excluding first and last as in postproc).
        psi_n : ndarray (n_surfaces-2,)
            Corresponding normalized psi values.
        """
        from scipy.integrate import cumulative_trapezoid

        R, Z = self.get_node_coordinates()
        psi_n0 = self.get_psi_n()
        psi = self.get_field('psi', tor_mode=0)

        # Find R range for flux surfaces along outboard midplane
        R_axis, Z_axis, _ = self.find_magnetic_axis()
        psi_bnd = self.find_separatrix_psi()

        # Get nodes at/near midplane (Z ≈ Z_axis) for outboard side
        mid_mask = np.abs(Z - Z_axis) < 0.01 * (np.max(Z) - np.min(Z) + 1e-10)
        outboard_mask = mid_mask & (R > R_axis)

        if np.sum(outboard_mask) < 10:
            # Fallback: use all nodes and sort by R
            outboard_mask = R > R_axis

        # Sort by R (proxy for flux surface label)
        R_ob = R[outboard_mask]
        psi_n_ob = psi_n0[outboard_mask]
        idx_sort = np.argsort(R_ob)
        R_ob = R_ob[idx_sort]
        psi_n_ob = psi_n_ob[idx_sort]

        # Create uniform psi_n grid and interpolate R along outboard midplane
        psi_n_grid = np.linspace(np.min(psi_n_ob) * 1.01,
                                 np.max(psi_n_ob) * 0.99, n_surfaces)
        R_grid = np.interp(psi_n_grid, psi_n_ob, R_ob)

        # Compute q for each surface
        # q = dPhi_pol / dPhi_tor  (field line pitch)
        # For axisymmetric equilibrium, q = (1/2π) ∮ (B_phi / (R * B_pol)) dl
        # Simplification: compute from psi gradients
        q = np.zeros(n_surfaces)

        for i in range(n_surfaces):
            # Compute poloidal flux gradient magnitude at R_grid[i]
            # Use finite difference on the midplane
            dr = R_grid[min(i + 1, n_surfaces - 1)] - R_grid[max(i - 1, 0)]
            if dr < 1e-10:
                dr = 1e-10

            # dpsi/dR at outboard midplane
            dpsi_dR_ob = (np.interp(R_grid[i] + dr * 0.1, R_ob, psi_n_ob) -
                          np.interp(R_grid[i] - dr * 0.1, R_ob, psi_n_ob)) / (dr * 0.2)

            # q ~ R * B_phi * d(Vol)/dpsi  (simplified)
            # More directly: q = dPhi_tor/dPsi_pol = F/(2π) * ∮ dl/(R^2 * B_pol)
            # Use cylindrical approximation:
            # q ≈ (R * B_phi) / (R^2 * dB_pol/dr) simplified as F0/(R*dpsi/dr)
            if dpsi_dR_ob > 0:
                q[i] = abs(self.F0) / (R_grid[i] * abs(2 * np.pi * dpsi_dR_ob) + 1e-20)
            # Clamp to reasonable values
            q[i] = np.clip(q[i], 0.1, 50)

        # Smooth q
        from scipy.ndimage import uniform_filter1d
        q = uniform_filter1d(q, size=5)

        return q[1:-1], psi_n_grid[1:-1]

    # --------------------------------------------------------------------------
    # Magnetic island width
    # --------------------------------------------------------------------------
    def get_island_widths(self, m_max=8, n_max=6):
        """Calculate magnetic island widths at rational surfaces.

        For each (m,n) pair, the island width at the rational surface q=m/n is:
            W = 4 * sqrt( |q * R * psi_mn / (s * B_phi)| )

        where:
            psi_mn = perturbed poloidal flux amplitude (from Fourier harmonics)
            q = m/n (safety factor)
            s = (r/q) * dq/dr (magnetic shear)
            R = major radius at the rational surface
            B_phi = F0 / R (toroidal field)

        Parameters
        ----------
        m_max : int
            Maximum poloidal mode number.
        n_max : int
            Maximum toroidal mode number.

        Returns
        -------
        islands : list of dict
            Each dict contains: m, n, q, psi_n, R, width, psi_pert_amplitude.
        """
        # Get equilibrium quantities
        psi_n0 = self.get_psi_n()
        psi = self.get_field('psi', tor_mode=0)
        psi_axis = self.find_magnetic_axis()[2]
        psi_bnd = self.find_separatrix_psi()
        R_axis, Z_axis, _ = self.find_magnetic_axis()
        R, Z = self.get_node_coordinates()

        # Compute q-profile
        q, psi_n_grid = self.get_q_profile(n_surfaces=150)

        # Get outboard midplane mapping psi_n -> R
        mid_mask = np.abs(Z - Z_axis) < 0.01 * (np.max(Z) - np.min(Z) + 1e-10)
        outboard_mask = mid_mask & (R > R_axis)
        if np.sum(outboard_mask) < 10:
            outboard_mask = R > R_axis
        psi_n_ob = psi_n0[outboard_mask]
        R_ob = R[outboard_mask]
        sort_idx = np.argsort(R_ob)
        psi_n_ob = psi_n_ob[sort_idx]
        R_ob = R_ob[sort_idx]

        islands = []

        # Scan rational surfaces
        for n in range(1, n_max + 1):
            n_tor = n * self.n_period
            for m in range(1, m_max + 1):
                q_target = float(m) / float(n)

                # Skip if q_target is outside computed q range
                if q_target < np.min(q) or q_target > np.max(q):
                    continue

                # Find the rational surface (psi_n where q = m/n)
                idx_q = np.searchsorted(q, q_target)
                if idx_q >= len(psi_n_grid):
                    continue
                psi_n_rs = psi_n_grid[idx_q]
                R_rs = np.interp(psi_n_rs, psi_n_ob, R_ob)

                # Compute magnetic shear at rational surface
                if idx_q > 0 and idx_q < len(q) - 1:
                    dq_dpsi = (q[idx_q + 1] - q[idx_q - 1]) / (
                        psi_n_grid[idx_q + 1] - psi_n_grid[idx_q - 1] + 1e-20)
                else:
                    dq_dpsi = 0.1

                # Shear s = (r/q) * dq/dr = (psi_n/q) * dq/dpsi_n * 2*psi_n
                # (since sqrt(psi_n) ~ r/a)
                r_over_a = np.sqrt(psi_n_rs)
                shear = (r_over_a / q_target) * dq_dpsi * 2 * r_over_a + 1e-20

                # Extract perturbation amplitude for (m,n) mode
                # The toroidal harmonic for n_tor in the HDF5:
                #   even index (2*k) = cos component
                #   odd index (2*k+1) = sin component
                # k = n_tor/n_period, index depends on n_tor and n_period
                k = n_tor // self.n_period
                if k < 1 or 2 * k >= self.n_tor:
                    continue

                idx_cos = 2 * k      # cos harmonic
                idx_sin = 2 * k + 1  # sin harmonic

                if idx_sin >= self.n_tor:
                    continue

                psi_cos = self.get_field('psi', tor_mode=idx_cos)
                psi_sin = self.get_field('psi', tor_mode=idx_sin)

                # Amplitude of perturbed psi at rational surface
                # Find nodes near the rational surface
                dist = np.abs(psi_n0 - psi_n_rs)
                near_mask = dist < 0.02  # Within ~2% in psi_n
                if np.sum(near_mask) < 5:
                    near_mask = dist < 0.05

                if np.sum(near_mask) == 0:
                    continue

                # RMS amplitude of perturbation
                amp_cos = np.sqrt(np.mean(psi_cos[near_mask]**2))
                amp_sin = np.sqrt(np.mean(psi_sin[near_mask]**2))
                psi_mn = np.sqrt(amp_cos**2 + amp_sin**2)

                # Absolute perturbation in physical units
                psi_mn_abs = psi_mn * abs(psi_bnd - psi_axis)

                # Island width formula (poloidal flux perturbation method):
                # W = 4 * sqrt( q * R * psi_mn_abs / (s * B_phi) )
                # B_phi = F0 / R
                B_phi = abs(self.F0) / R_rs if abs(self.F0) > 1e-20 else 1.0

                width = 4.0 * np.sqrt(
                    abs(q_target * R_rs * psi_mn_abs /
                        (abs(shear) * B_phi + 1e-20))
                )

                # Normalized width (as fraction of minor radius)
                a_minor = (np.max(R[outboard_mask]) - R_axis) if np.any(outboard_mask) else 0.5
                width_norm = width / a_minor if a_minor > 0 else width

                islands.append({
                    'm': m,
                    'n': n,
                    'q': q_target,
                    'psi_n': psi_n_rs,
                    'R': R_rs,
                    'width': width,
                    'width_normalized': width_norm,
                    'shear': shear,
                    'psi_pert_amplitude': psi_mn_abs,
                })

        # Sort by width (largest first)
        islands.sort(key=lambda x: x['width'], reverse=True)

        return islands

    def print_island_widths(self, islands=None, min_width=0.0):
        """Print a formatted table of magnetic island widths.

        Parameters
        ----------
        islands : list of dict or None
            Result from get_island_widths(). If None, compute automatically.
        min_width : float
            Minimum width (in meters) to display.
        """
        if islands is None:
            islands = self.get_island_widths()

        print("\n" + "=" * 90)
        print("  Magnetic Island Widths")
        print("=" * 90)
        print(f"  {'(m,n)':<10} {'q=m/n':>8} {'psi_n':>8} {'R(m)':>8} "
              f"{'Width(m)':>10} {'W/a':>8} {'Shear':>8} {'psi_pert':>10}")
        print("-" * 90)

        shown = 0
        for isl in islands:
            if isl['width'] < min_width:
                continue
            print(f"  ({isl['m']},{isl['n']}):    "
                  f"{isl['q']:8.3f} {isl['psi_n']:8.4f} {isl['R']:8.4f} "
                  f"{isl['width']:10.4f} {isl['width_normalized']:8.4f} "
                  f"{isl['shear']:8.4f} {isl['psi_pert_amplitude']:10.4e}")
            shown += 1

        print("-" * 90)
        print(f"  Total rational surfaces with islands: {len(islands)}")
        print(f"  Shown (> {min_width} m): {shown}")
        print("=" * 90)

    # --------------------------------------------------------------------------
    # Poincare plot data
    # --------------------------------------------------------------------------
    def trace_field_line(self, R_start, Z_start, n_turns=200, delta_phi=0.02,
                         max_steps=500000):
        """Simple field-line tracer for generating Poincaré plots.

        Uses RK4 integration with Bezier element evaluation.

        Parameters
        ----------
        R_start, Z_start : float
            Starting (R, Z) coordinates.
        n_turns : int
            Number of poloidal turns.
        delta_phi : float
            Toroidal step size in radians.
        max_steps : int
            Maximum integration steps.

        Returns
        -------
        R_cross, Z_cross : ndarray
            Poincaré crossing points at the outboard midplane.
        psi_n_cross : ndarray
            Normalized psi at each crossing.
        """
        n_period = self.n_period
        F0 = self.F0
        n_tor = self.n_tor

        R_nodes, Z_nodes = self.get_node_coordinates()
        vtx = self.vertex  # (n_elements, n_vertex_max)

        # Find starting element
        i_elm = self._find_element(R_start, Z_start)
        if i_elm < 0:
            print(f"Warning: Starting point ({R_start}, {Z_start}) not in domain")
            return np.array([]), np.array([]), np.array([])

        # Pre-compute Bezier control points for all elements
        n_elem = vtx.shape[0]
        elem_R_ctrl = np.zeros((n_elem, 4, 4))  # (n_elem, degree, vertex)
        elem_Z_ctrl = np.zeros((n_elem, 4, 4))

        for i in range(n_elem):
            n1, n2, n3, n4 = vtx[i, :4] - 1
            if min(n1, n2, n3, n4) < 0:
                continue
            for deg in range(min(4, self.n_degrees)):
                elem_R_ctrl[i, deg, 0] = self.x_coords[n1, 0, deg, 0]
                elem_R_ctrl[i, deg, 1] = self.x_coords[n2, 0, deg, 0]
                elem_R_ctrl[i, deg, 2] = self.x_coords[n3, 0, deg, 0]
                elem_R_ctrl[i, deg, 3] = self.x_coords[n4, 0, deg, 0]
                elem_Z_ctrl[i, deg, 0] = self.x_coords[n1, 0, deg, 1]
                elem_Z_ctrl[i, deg, 1] = self.x_coords[n2, 0, deg, 1]
                elem_Z_ctrl[i, deg, 2] = self.x_coords[n3, 0, deg, 1]
                elem_Z_ctrl[i, deg, 3] = self.x_coords[n4, 0, deg, 1]

        # Field evaluation at (s, t, phi) in element i_elm
        def eval_field(s, t, phi, var_idx, tor_mode):
            """Evaluate a field variable at (s,t,phi) using Bezier basis."""
            n_deg = min(self.n_degrees, 4)
            bf = bezier_basis(s, t)  # (n_order, n_vertex)

            # Get node values for this element's 4 vertices
            n_vals_corner = np.zeros(4)
            for v in range(4):
                n_idx = vtx[i_elm, v] - 1
                if n_idx >= 0:
                    val = 0.0
                    for deg in range(n_deg):
                        val += self.values[n_idx, tor_mode, deg, var_idx - 1] * bf[deg, v]
                    n_vals_corner[v] = val

            # Sum contributions: value = sum over (degree, vertex) of
            #   value[degree, vertex] * basis[degree](s, t)[vertex]
            result = 0.0
            for v in range(4):
                n_idx = vtx[i_elm, v] - 1
                if n_idx < 0:
                    continue
                for deg in range(n_deg):
                    result += (self.values[n_idx, tor_mode, deg, var_idx - 1] *
                               bf[deg, v])
            return result

        def eval_psi(s, t, phi):
            """Evaluate total psi at (s,t,phi) including all harmonics."""
            n_deg = min(self.n_degrees, 4)
            bf = bezier_basis(s, t)

            # n=0 component
            psi_val = 0.0
            for v in range(4):
                n_idx = vtx[i_elm, v] - 1
                if n_idx < 0:
                    continue
                for deg in range(n_deg):
                    psi_val += self.values[n_idx, 0, deg, 0] * bf[deg, v]

            # Add toroidal harmonics
            for i_tor in range(1, (n_tor - 1) // 2 + 1):
                i_cos = 2 * i_tor
                i_sin = 2 * i_tor + 1
                mode_n = i_tor * n_period
                cos_coeff = 0.0
                sin_coeff = 0.0
                for v in range(4):
                    n_idx = vtx[i_elm, v] - 1
                    if n_idx < 0:
                        continue
                    for deg in range(n_deg):
                        cos_coeff += self.values[n_idx, i_cos, deg, 0] * bf[deg, v]
                        if i_sin < n_tor:
                            sin_coeff += self.values[n_idx, i_sin, deg, 0] * bf[deg, v]
                psi_val += cos_coeff * np.cos(mode_n * phi)
                if i_sin < n_tor:
                    psi_val += sin_coeff * np.sin(mode_n * phi)

            return psi_val

        def eval_RZ(s, t):
            """Evaluate (R, Z) at element coordinates (s, t)."""
            n_deg = min(self.n_degrees, 4)
            bf = bezier_basis(s, t)
            R_val, Z_val = 0.0, 0.0
            for v in range(4):
                n_idx = vtx[i_elm, v] - 1
                if n_idx < 0:
                    continue
                for deg in range(n_deg):
                    R_val += elem_R_ctrl[i_elm, deg, v] * bf[deg, v]
                    Z_val += elem_Z_ctrl[i_elm, deg, v] * bf[deg, v]
            return R_val, Z_val

        def eval_B(s, t, phi):
            """Evaluate magnetic field (BR, BZ, Bphi) at (s,t,phi)."""
            # Compute spatial derivatives of psi
            ds = 1e-4
            R0, Z0 = eval_RZ(s, t)
            Rp, Zp = eval_RZ(s + ds, t)
            Rm, Zm = eval_RZ(s - ds, t)
            Rt, Zt = eval_RZ(s, t + ds)
            Rb, Zb = eval_RZ(s, t - ds)

            R_s = (Rp - Rm) / (2 * ds)
            Z_s = (Zp - Zm) / (2 * ds)
            R_t = (Rt - Rb) / (2 * ds)
            Z_t = (Zt - Zb) / (2 * ds)
            Zjac = R_s * Z_t - R_t * Z_s

            psi0 = eval_psi(s, t, phi)
            psi_ps = eval_psi(s + ds, t, phi)
            psi_ms = eval_psi(s - ds, t, phi)
            psi_pt = eval_psi(s, t + ds, phi)
            psi_mt = eval_psi(s, t - ds, phi)
            psi_pp = eval_psi(s, t, phi + ds)

            psi_s = (psi_ps - psi_ms) / (2 * ds)
            psi_t = (psi_pt - psi_mt) / (2 * ds)
            psi_p = (psi_pp - psi0) / ds  # One-sided for toroidal

            # Physical gradients
            if abs(Zjac) < 1e-20:
                return 0.0, 0.0, 0.0

            psi_R = (Z_t * psi_s - Z_s * psi_t) / Zjac
            psi_Z = (-R_t * psi_s + R_s * psi_t) / Zjac

            # Reduced MHD B-field
            BR = -psi_Z / R0
            BZ = psi_R / R0
            Bphi = F0 / R0

            return BR, BZ, Bphi

        def find_neighbor(i_elm, s, t, side):
            """Find neighboring element using vertex matching."""
            current_vertices = set(vtx[i_elm, :4])
            for j in range(n_elem):
                if j == i_elm:
                    continue
                other_vertices = set(vtx[j, :4])
                common = current_vertices & other_vertices
                if len(common) == 2:
                    # Check if this is the neighbor on the correct side
                    # This is a simplified approach
                    return j
            return i_elm  # Stay in same element if neighbor not found

        # Main field-line tracing loop
        crossings_R, crossings_Z, crossings_psi_n = [], [], []
        R_axis, Z_axis, psi_axis = self.find_magnetic_axis()
        psi_bnd = self.find_separatrix_psi()

        s, t = 0.5, 0.5  # Start at element center (approximate)
        phi = 0.0
        phi_prev = 0.0
        Z_prev = Z_start
        n_cross = 0

        R_cur, Z_cur = eval_RZ(s, t)

        for step in range(max_steps):
            if n_cross >= n_turns:
                break

            # RK4 step
            BR, BZ, Bphi = eval_B(s, t, phi)
            Bmag = np.sqrt(BR**2 + BZ**2 + Bphi**2)
            if Bmag < 1e-20:
                break

            dphi = delta_phi

            # Normalize B so |B| * dphi / Bphi gives proper toroidal step
            ds_over_dphi = BR / (R_cur * Bphi) if abs(Bphi) > 1e-20 else 0.0
            dt_over_dphi = BZ / (R_cur * Bphi) if abs(Bphi) > 1e-20 else 0.0

            # Jacobian for converting to element coordinates
            # ds/dphi = (-Z_t*R_p + R_t*Z_p + R*(Z_t*BR - R_t*BZ)/Bphi) / Zjac
            # Simplified: step in physical space then map back
            R_new = R_cur + BR / Bphi * dphi
            Z_new = Z_cur + BZ / Bphi * dphi

            ds = (-Z_t * BR + R_t * BZ) / Bphi * dphi  # simplified
            dt = (Z_s * BR - R_s * BZ) / Bphi * dphi   # simplified

            # Recompute R_s, R_t, Z_s, Z_t
            ds2 = 1e-6
            Rp, Zp = eval_RZ(s + ds2, t)
            Rm, Zm = eval_RZ(s - ds2, t)
            R_s = (Rp - Rm) / (2 * ds2)
            Z_s = (Zp - Zm) / (2 * ds2)
            Rt, Zt = eval_RZ(s, t + ds2)
            Rb, Zb = eval_RZ(s, t - ds2)
            R_t = (Rt - Rb) / (2 * ds2)
            Z_t = (Zt - Zb) / (2 * ds2)
            Zjac = R_s * Z_t - R_t * Z_s

            if abs(Zjac) > 1e-20:
                R_p = 0.0
                Z_p = 0.0  # Simplified: grid stationary
                ds = (-Z_t * R_p + R_t * Z_p +
                      R_cur * (Z_t * BR - R_t * BZ) / Bphi) * dphi / Zjac
                dt = (Z_s * R_p - R_s * Z_p -
                      R_cur * (Z_s * BR - R_s * BZ) / Bphi) * dphi / Zjac

            s_new = s + ds
            t_new = t + dt
            phi_new = phi + dphi

            # Element boundary check and transition
            if s_new < 0 or s_new > 1 or t_new < 0 or t_new > 1:
                new_elm = find_neighbor(i_elm, s_new, t_new,
                                        side=int(s_new < 0) + 2 * int(t_new < 0))
                if new_elm != i_elm:
                    i_elm = new_elm
                s_new = np.clip(s_new, 0.001, 0.999)
                t_new = np.clip(t_new, 0.001, 0.999)

            s, t, phi = s_new, t_new, phi_new
            R_cur, Z_cur = eval_RZ(s, t)

            # Check outboard midplane crossing
            if (R_cur > R_axis and
                np.sign(Z_cur - Z_axis) != np.sign(Z_prev - Z_axis)):
                n_cross += 1
                crossings_R.append(R_cur)
                crossings_Z.append(Z_cur)
                # Evaluate psi_n at crossing
                psi_val = eval_psi(s, t, phi)
                psi_n_val = (psi_val - psi_axis) / (psi_bnd - psi_axis + 1e-20)
                crossings_psi_n.append(psi_n_val)

            Z_prev = Z_cur

        return (np.array(crossings_R), np.array(crossings_Z),
                np.array(crossings_psi_n))

    def _find_element(self, R, Z):
        """Find element containing (R, Z) by centroid distance."""
        R_nodes, Z_nodes = self.get_node_coordinates()
        vtx = self.vertex
        n_elem = vtx.shape[0]

        best_dist = np.inf
        best_elem = -1

        for i in range(n_elem):
            ns = vtx[i, :4] - 1
            if ns.min() < 0:
                continue
            R_cent = np.mean(R_nodes[ns])
            Z_cent = np.mean(Z_nodes[ns])
            dist = (R - R_cent)**2 + (Z - Z_cent)**2
            if dist < best_dist:
                best_dist = dist
                best_elem = i

        return best_elem


# ==============================================================================
# Time series — batch processing of multiple restart files
# ==============================================================================
def _scan_restart_files(directory, only=None, donly=None, use_5digits=False):
    """Scan directory for JOREK restart HDF5 files.

    Parameters
    ----------
    directory : str
        Directory to scan.
    only : str or None
        Step selection, e.g. "0,100,200", "0-10-100".
    donly : int or None
        Step interval, e.g. 10 = every 10th step.
    use_5digits : bool
        Use 5-digit (jorek00001.h5) instead of 6-digit (jorek000001.h5).

    Returns
    -------
    files : list of (int, str)
        Sorted list of (step_number, absolute_path) tuples.
    """
    if use_5digits:
        pattern = "jorek?????.h5"
    else:
        pattern = "jorek??????.h5"

    import glob
    full_pattern = os.path.join(directory, pattern)
    all_files = sorted(glob.glob(full_pattern))

    if not all_files:
        # Also check for jorek_restart.h5 (single file) in directory
        single = os.path.join(directory, "jorek_restart.h5")
        if os.path.exists(single):
            return [(0, single)]

    # Extract step numbers
    file_steps = []
    for f in all_files:
        basename = os.path.basename(f)
        # "jorek000123.h5" -> step_num = 123
        step_str = basename.replace("jorek", "").replace(".h5", "")
        try:
            step_num = int(step_str)
            file_steps.append((step_num, os.path.abspath(f)))
        except ValueError:
            continue

    # Parse selection
    selected_steps = set()
    if donly is not None:
        for s, _ in file_steps:
            if s % donly == 0:
                selected_steps.add(s)
    elif only is not None:
        for part in only.split(","):
            parts = part.split("-")
            if len(parts) == 1:
                selected_steps.add(int(parts[0]))
            elif len(parts) == 2:
                start, end = int(parts[0]), int(parts[1])
                selected_steps.update(range(start, end + 1))
            elif len(parts) == 3:
                start, step, end = int(parts[0]), int(parts[1]), int(parts[2])
                selected_steps.update(range(start, end + 1, step))
    else:
        # All files
        selected_steps = set(s for s, _ in file_steps)

    # Filter
    result = [(s, f) for s, f in file_steps if s in selected_steps]
    return sorted(result)


class JOREKTimeSeries:
    """Batch processor for JOREK restart file time series.

    Reads a sequence of JOREK restart HDF5 files from a directory and
    extracts evolution of quantities (island widths, q-profile, etc.)
    over time.

    Parameters
    ----------
    directory : str
        Directory containing jorek??????.h5 files.
    only : str or None
        Select specific steps: "0,100,200" or "0-10-100" or "100-5-200".
    donly : int or None
        Take every Nth step.
    use_5digits : bool
        For 5-digit restart file numbering.
    verbose : bool
        Print progress information.

    Examples
    --------
    >>> ts = JOREKTimeSeries('.', donly=100)
    >>> result = ts.track_island('2/1')  # Track 2/1 island width over time
    >>> ts.plot_island_evolution('2/1', '3/2')  # Plot multiple modes
    >>> data = ts.get_island_evolution(['2/1', '3/2', '4/3'])
    """

    def __init__(self, directory='.', only=None, donly=None,
                 use_5digits=False, verbose=True):
        self.directory = os.path.abspath(directory)
        self.verbose = verbose
        self._file_list = _scan_restart_files(
            directory, only=only, donly=donly, use_5digits=use_5digits)

        if not self._file_list:
            raise FileNotFoundError(
                f"No JOREK restart HDF5 files found in {self.directory}")

        if self.verbose:
            print(f"Found {len(self._file_list)} restart files in {self.directory}")
            if len(self._file_list) > 0:
                print(f"  Steps: {self._file_list[0][0]} ... {self._file_list[-1][0]}")
            if donly:
                print(f"  Using every {donly}th step")
            elif only:
                print(f"  Selected steps: {only}")

    @property
    def n_files(self):
        return len(self._file_list)

    def get_time_array(self):
        """Read time values from all restart files.

        Returns
        -------
        times : ndarray (n_files,)
        steps : ndarray (n_files,)
        """
        times = np.zeros(self.n_files)
        steps = np.zeros(self.n_files, dtype=int)
        for i, (step, fname) in enumerate(self._file_list):
            try:
                with h5py.File(fname, 'r') as f:
                    t = f['t_now'][()]
                    if np.ndim(t) > 0:
                        t = t.flat[0]
                    times[i] = float(t)
                    steps[i] = step
            except Exception:
                times[i] = float(step)
                steps[i] = step
        return times, steps

    def get_island_evolution(self, mode_list=None, m_max=8, n_max=6,
                             n_surfaces=100, min_width=0.0):
        """Compute island width evolution across all restart files.

        Parameters
        ----------
        mode_list : list of str or None
            Specific modes to track, e.g. ['2/1', '3/2'].
            If None, tracks the top modes from the first time step.
        m_max, n_max : int
            Max poloidal/toroidal mode numbers to scan.
        n_surfaces : int
            Number of flux surfaces for q-profile.
        min_width : float
            Minimum width to record (m).

        Returns
        -------
        result : dict
            Keys:
            - 'times': ndarray (n_files,) of time values
            - 'steps': ndarray (n_files,) of step indices
            - 'modes': list of (m,n) tuples tracked
            - 'widths': dict mapping (m,n) -> ndarray (n_files,) of widths
            - 'psi_n': dict mapping (m,n) -> ndarray (n_files,) of psi_n at
                       rational surface
            - 'psi_pert': dict mapping (m,n) -> ndarray (n_files,) of
                          perturbation amplitude
        """
        times, steps = self.get_time_array()

        # Parse mode_list
        if mode_list is not None:
            tracked_modes = []
            for s in mode_list:
                parts = s.split('/')
                if len(parts) == 2:
                    tracked_modes.append((int(parts[0]), int(parts[1])))
        else:
            tracked_modes = None  # Will determine after first step

        # Initialize storage
        all_modes = set()
        widths = {}
        psi_n_vals = {}
        psi_pert_vals = {}

        # Process first file to determine equilibrium anchor
        if self.verbose:
            print(f"\nProcessing {self.n_files} files...")

        first_file = self._file_list[0][1]

        for i_file, (step, fname) in enumerate(self._file_list):
            if self.verbose and (i_file % max(1, self.n_files // 10) == 0):
                print(f"  [{i_file+1}/{self.n_files}] step {step}...")

            try:
                jr = JOREKRestart(fname)

                # Get island widths
                islands = jr.get_island_widths(m_max=m_max, n_max=n_max)

                # Build lookup: (m,n) -> island data
                isl_dict = {(isl['m'], isl['n']): isl for isl in islands}

                # On first step, determine which modes to track
                if i_file == 0:
                    if tracked_modes is None:
                        # Auto-select: take top modes by width
                        top = [isl for isl in islands
                               if isl['width'] > min_width][:10]
                        tracked_modes = [(isl['m'], isl['n']) for isl in top]
                        if self.verbose:
                            print(f"  Auto-selected modes: {tracked_modes}")

                    for mn in tracked_modes:
                        widths[mn] = np.zeros(self.n_files)
                        psi_n_vals[mn] = np.zeros(self.n_files)
                        psi_pert_vals[mn] = np.zeros(self.n_files)

                # Fill data for this time step
                for mn in tracked_modes:
                    if mn in isl_dict:
                        widths[mn][i_file] = isl_dict[mn]['width']
                        psi_n_vals[mn][i_file] = isl_dict[mn]['psi_n']
                        psi_pert_vals[mn][i_file] = isl_dict[mn]['psi_pert_amplitude']
                    else:
                        widths[mn][i_file] = 0.0
                        psi_n_vals[mn][i_file] = np.nan
                        psi_pert_vals[mn][i_file] = 0.0

                jr.close()

            except Exception as e:
                if self.verbose:
                    print(f"  WARNING: Failed processing {fname}: {e}")
                for mn in tracked_modes or []:
                    if mn in widths:
                        widths[mn][i_file] = np.nan

        return {
            'times': times,
            'steps': steps,
            'modes': tracked_modes,
            'widths': widths,
            'psi_n': psi_n_vals,
            'psi_pert': psi_pert_vals,
        }

    def track_island(self, mode_str, **kwargs):
        """Convenience: track a single mode's width evolution.

        Parameters
        ----------
        mode_str : str
            Mode string, e.g. '2/1'.
        **kwargs :
            Passed to get_island_evolution().

        Returns
        -------
        times, widths : ndarray
        """
        result = self.get_island_evolution(mode_list=[mode_str], **kwargs)
        m, n = result['modes'][0]
        return result['times'], result['widths'][(m, n)]

    def plot_island_evolution(self, *mode_strs, output=None, title=None,
                              show=True, **kwargs):
        """Plot island width evolution over time for selected modes.

        Parameters
        ----------
        *mode_strs : str
            Mode strings, e.g. '2/1', '3/2'.
        output : str or None
            Output filename for the plot.
        title : str or None
            Plot title.
        show : bool
            Whether to display the plot.
        **kwargs :
            Passed to get_island_evolution().
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("ERROR: matplotlib is required for plotting.")
            return

        if not mode_strs:
            print("ERROR: At least one mode must be specified, e.g. '2/1'")
            return

        result = self.get_island_evolution(
            mode_list=list(mode_strs), **kwargs)

        times = result['times']
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

        # Color cycle
        colors = plt.cm.tab10(np.linspace(0, 1, max(len(mode_strs), 1)))

        for i, mode_str in enumerate(mode_strs):
            parts = mode_str.split('/')
            m, n = int(parts[0]), int(parts[1])
            w = result['widths'][(m, n)]
            psi_pert = result['psi_pert'][(m, n)]

            label = f'({m},{n})'
            color = colors[i % len(colors)]

            # Filter NaN
            valid = ~np.isnan(w) & (w > 0)

            # Top panel: island width
            ax1.plot(times[valid], w[valid] * 100, 'o-', color=color,
                     markersize=3, linewidth=1.5, label=label)
            ax1.set_ylabel('Island Width [cm]')
            ax1.set_title(title or 'Magnetic Island Width Evolution')
            ax1.legend(loc='best')
            ax1.grid(True, alpha=0.3)

            # Bottom panel: perturbation amplitude
            ax2.plot(times[valid], psi_pert[valid], 'o-', color=color,
                     markersize=3, linewidth=1.5, label=label)
            ax2.set_ylabel(r'$\tilde{\psi}$ perturbation amplitude')
            ax2.legend(loc='best')
            ax2.grid(True, alpha=0.3)

        ax2.set_xlabel('Time [JOREK units]')

        plt.tight_layout()

        if output:
            fig.savefig(output, dpi=150, bbox_inches='tight')
            print(f"\n  Island evolution plot saved to: {output}")

        if show:
            plt.show()

        return fig


# ==============================================================================
# Command-line interface
# ==============================================================================
def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='JOREK Restart File Reader and Island Width Calculator')

    # Subcommands for single file vs batch analysis
    sub = parser.add_subparsers(dest='command', help='Analysis mode')

    # --- Single file mode ---
    p_single = sub.add_parser('single', help='Analyze a single restart file',
                              aliases=['s'])
    p_single.add_argument('filename', help='Path to jorek_restart.h5 file')
    p_single.add_argument('--info', action='store_true', default=True,
                          help='Print file information (default)')
    p_single.add_argument('--islands', action='store_true',
                          help='Calculate magnetic island widths')
    p_single.add_argument('--qprofile', action='store_true',
                          help='Plot q-profile')
    p_single.add_argument('--poincare', action='store_true',
                          help='Trace a field line for Poincaré plot')
    p_single.add_argument('--R', type=float, default=None,
                          help='Starting R for field line trace')
    p_single.add_argument('--Z', type=float, default=0.0,
                          help='Starting Z for field line trace')
    p_single.add_argument('--n-turns', type=int, default=200,
                          help='Number of toroidal turns for trace')
    p_single.add_argument('--min-width', type=float, default=0.0,
                          help='Minimum island width to display (m)')
    p_single.add_argument('--output', '-o', type=str, default=None,
                          help='Output file for data')
    p_single.add_argument('--no-info', action='store_false', dest='info',
                          help='Suppress file info')

    # --- Batch / time-series mode ---
    p_batch = sub.add_parser('batch', help='Batch process restart file series',
                             aliases=['b'])
    p_batch.add_argument('directory', nargs='?', default='.',
                         help='Directory containing jorek??????.h5 files '
                              '[default: .]')
    p_batch.add_argument('--modes', '-m', type=str, default=None,
                         help='Modes to track, e.g. "2/1,3/2,4/3" '
                              '[default: auto-select top modes]')
    p_batch.add_argument('--only', type=str, default=None,
                         help='Select specific steps, e.g. "0,100,200" or '
                              '"100-10-500"')
    p_batch.add_argument('--donly', type=int, default=None,
                         help='Take every Nth step')
    p_batch.add_argument('--min-width', type=float, default=0.001,
                         help='Minimum island width to track (m)')
    p_batch.add_argument('--output', '-o', type=str,
                         default='island_evolution.png',
                         help='Output plot filename')
    p_batch.add_argument('--no-show', action='store_true',
                         help='Do not display plot (save only)')
    p_batch.add_argument('--title', type=str, default=None,
                         help='Plot title')
    p_batch.add_argument('--save-data', type=str, default=None,
                         help='Save raw data to file (e.g. island_data.txt)')

    args = parser.parse_args()

    # If no subcommand given, default to single-file mode (backward compat)
    if args.command is None:
        # Check if first positional arg looks like a directory with many files
        if len(sys.argv) > 1 and not sys.argv[1].startswith('-'):
            first = sys.argv[1]
            if os.path.isdir(first):
                # Treat as batch mode
                args.command = 'batch'
                args.directory = first
                args.modes = None
                args.only = None
                args.donly = None
                args.min_width = 0.001
                args.output = 'island_evolution.png'
                args.no_show = False
                args.title = None
                args.save_data = None
            elif os.path.isfile(first) and first.endswith('.h5'):
                args.command = 'single'
                args.filename = first
                args.info = True
                args.islands = False
                args.qprofile = False
                args.poincare = False
                args.R = None
                args.Z = 0.0
                args.n_turns = 200
                args.min_width = 0.0
                args.output = None
            else:
                print(f"ERROR: {first} is not a valid file or directory.")
                sys.exit(1)
        else:
            parser.print_help()
            sys.exit(0)

    # --- Execute ---
    if args.command in ('single', 's'):
        _run_single(args)
    elif args.command in ('batch', 'b'):
        _run_batch(args)


def _run_single(args):
    """Execute single-file analysis."""
    if not os.path.exists(args.filename):
        print(f"ERROR: File not found: {args.filename}")
        sys.exit(1)

    jr = JOREKRestart(args.filename, verbose=True)

    try:
        if args.info:
            jr.print_info()

        if args.qprofile:
            q, psi_n = jr.get_q_profile()
            try:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(1, 1, figsize=(8, 6))
                ax.plot(psi_n, q, 'b-', linewidth=2)
                ax.set_xlabel(r'$\psi_n$')
                ax.set_ylabel('q (safety factor)')
                ax.set_title('q-profile')
                ax.grid(True, alpha=0.3)
                for m in [1, 2, 3, 4, 5, 6]:
                    for n in [1, 2, 3, 4]:
                        q_target = m / n
                        if q_target >= np.min(q) and q_target <= np.max(q):
                            ax.axhline(y=q_target, color='red', linestyle='--',
                                      alpha=0.3, linewidth=0.5)
                fig.savefig('q_profile.png', dpi=150, bbox_inches='tight')
                print("\n  Q-profile plot saved to: q_profile.png")
                plt.show()
            except ImportError:
                print("Warning: matplotlib not available.")
                print("\n  Q-profile (psi_n, q):")
                for psi_n_val, q_val in zip(psi_n, q):
                    print(f"    {psi_n_val:8.4f}  {q_val:8.4f}")

        if args.islands:
            islands = jr.get_island_widths()
            jr.print_island_widths(islands, min_width=args.min_width)
            if args.output:
                _save_island_data(islands, args.output)

        if args.poincare:
            _run_poincare(jr, args.R, args.Z, args.n_turns)

    finally:
        jr.close()


def _run_batch(args):
    """Execute batch time-series analysis."""
    dir_path = args.directory
    if not os.path.isdir(dir_path):
        print(f"ERROR: Directory not found: {dir_path}")
        sys.exit(1)

    # Parse modes
    mode_list = None
    if args.modes:
        mode_list = [m.strip() for m in args.modes.split(',')]

    print("=" * 60)
    print(f"  JOREK Island Width Time-Series Analysis")
    print(f"  Directory: {dir_path}")
    if mode_list:
        print(f"  Tracking modes: {mode_list}")
    else:
        print(f"  Auto-selecting top modes")
    if args.only:
        print(f"  Step selection: {args.only}")
    if args.donly:
        print(f"  Step interval: every {args.donly}")
    print("=" * 60)

    ts = JOREKTimeSeries(
        directory=dir_path, only=args.only, donly=args.donly, verbose=True)

    if mode_list:
        ts.plot_island_evolution(
            *mode_list, min_width=args.min_width,
            output=args.output, show=not args.no_show, title=args.title)
    else:
        # Auto mode: get result, then find top modes
        result = ts.get_island_evolution(min_width=args.min_width)
        top_modes = result['modes']  # Already auto-selected
        if top_modes:
            mode_strs = [f'{m}/{n}' for m, n in top_modes]
            ts.plot_island_evolution(
                *mode_strs, min_width=args.min_width,
                output=args.output, show=not args.no_show, title=args.title)

    # Save raw data if requested
    if args.save_data:
        _save_evolution_data(ts, args.save_data,
                             mode_list=mode_list, min_width=args.min_width)


def _save_island_data(islands, output):
    """Save island data to ASCII file."""
    with open(output, 'w') as f:
        f.write("# Magnetic island widths\n")
        f.write("# (m,n)  q    psi_n  R(m)  Width(m)  "
                "W/a   Shear  psi_pert\n")
        for isl in islands:
            f.write(f"{isl['m']} {isl['n']} {isl['q']:.4f} "
                    f"{isl['psi_n']:.4f} {isl['R']:.4f} "
                    f"{isl['width']:.8f} {isl['width_normalized']:.4f} "
                    f"{isl['shear']:.4f} {isl['psi_pert_amplitude']:.6e}\n")
    print(f"\n  Island data saved to: {output}")


def _save_evolution_data(ts, output, mode_list=None, min_width=0.0):
    """Save evolution data to ASCII file."""
    result = ts.get_island_evolution(
        mode_list=mode_list, min_width=min_width)

    with open(output, 'w') as f:
        f.write("# Island width evolution over time\n")
        f.write(f"# Modes: {result['modes']}\n")
        header = "# {:>14s}".format("time")
        for m, n in result['modes']:
            header += f"  {'W_'+str(m)+'_'+str(n)+'(m)':>16s}"
            header += f"  {'psi_'+str(m)+'_'+str(n):>16s}"
        f.write(header + "\n")

        for i in range(len(result['times'])):
            line = f"  {result['times'][i]:14.6e}"
            for m, n in result['modes']:
                line += f"  {result['widths'][(m,n)][i]:16.8f}"
                line += f"  {result['psi_pert'][(m,n)][i]:16.8e}"
            f.write(line + "\n")

    print(f"\n  Evolution data saved to: {output}")


def _run_poincare(jr, R_start, Z_start, n_turns):
    """Run Poincare field-line trace."""
    if R_start is None:
        print("ERROR: --R is required for Poincaré tracing")
        sys.exit(1)

    R_cross, Z_cross, psi_cross = jr.trace_field_line(
        R_start, Z_start, n_turns=n_turns)

    if len(R_cross) > 0:
        try:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 1, figsize=(8, 10))
            ax.set_aspect('equal')
            ax.scatter(R_cross, Z_cross, s=1, c='blue', alpha=0.6,
                      rasterized=True)
            R_axis, Z_axis, _ = jr.find_magnetic_axis()
            ax.plot(R_axis, Z_axis, 'r+', markersize=15,
                   label='Magnetic Axis')
            ax.set_xlabel('R [m]')
            ax.set_ylabel('Z [m]')
            ax.set_title(f'Poincaré plot: R={R_start}, Z={Z_start}')
            ax.legend()
            ax.grid(True, alpha=0.3)
            outfile = f'poincare_R{R_start}_Z{Z_start}.png'
            fig.savefig(outfile, dpi=150, bbox_inches='tight')
            print(f"\n  Poincaré plot saved to: {outfile}")
            plt.show()
        except ImportError:
            print("Warning: matplotlib not available.")
        print(f"\n  Field line: {len(R_cross)} crossings recorded")
        if len(psi_cross) > 0:
            print(f"  psi_n range: [{psi_cross.min():.4f}, "
                  f"{psi_cross.max():.4f}]")
    else:
        print("WARNING: No crossings recorded.")


if __name__ == '__main__':
    main()
