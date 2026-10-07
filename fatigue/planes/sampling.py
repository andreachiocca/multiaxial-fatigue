# fatigue/planes/sampling.py
from __future__ import annotations
import numpy as np
from typing import Union, Optional


def dcmesh_xyz(a1, n1, a2, n2=None):
    """
    Generate points on a unit-sphere patch (direction cosines), translated from the MATLAB dcmesh().

    Parameters
    ----------
    a1 : float or array-like of length 2
        Angular span of interest in the XOZ plane (radians).
        If scalar: interpreted as [0, a1].
        If length-2: interpreted as [a1_min, a1_max].
        Must satisfy 0 <= a1 <= pi.

    n1 : int
        Number of subdivisions of the a1 range.

    a2 : float or array-like of length 2
        Angular span of interest in the XOY plane (radians).
        If scalar: interpreted as [0, a2].
        If length-2: interpreted as [a2_min, a2_max].
        Must satisfy 0 <= a2 <= 2*pi.

    n2 : int or None, optional
        If given, number of subdivisions of the a2 range used to set the target angular step da2.
        If None, da2 is set equal to da1 (as in the MATLAB code).

    Returns
    -------
    xyz : (M, 3) numpy.ndarray
        Points on the unit sphere (radius = 1). Each row is [x, y, z].
        Since radius = 1, these are direction cosines.

    Notes
    -----
    - This function intentionally does NOT compute triangulation (Delaunay) and does not plot anything.
    - Like the MATLAB version, if a2 spans a full circle (e.g. 0..2*pi), each latitude ring includes
      both endpoints, which yields a duplicate point at phi=0 and phi=2*pi. If you want to remove
      duplicates, you can post-process xyz (see note at end of this message).
    """

    # -----------------------------
    # Helper: parse scalar vs [min,max]
    # -----------------------------
    def _parse_angle_range(a, name, amax):
        arr = np.asarray(a, dtype=float)

        if arr.ndim == 0:
            lo, hi = 0.0, float(arr)
        else:
            arr = arr.ravel()
            if arr.size == 1:
                lo, hi = 0.0, float(arr[0])
            elif arr.size == 2:
                lo, hi = float(arr[0]), float(arr[1])
            else:
                raise ValueError(f"{name} must be a scalar or a length-2 sequence.")

        # MATLAB code checks elementwise bounds; here we check both endpoints
        if lo < 0 or hi < 0 or lo > amax or hi > amax:
            raise ValueError(f"Range ERROR: 0 <= {name} <= {amax}")

        if hi < lo:
            raise ValueError(f"Invalid range for {name}: max must be >= min.")

        return lo, hi

    # -----------------------------
    # Basic input validation
    # -----------------------------
    n1 = int(n1)
    if n1 <= 0:
        raise ValueError("n1 must be a positive integer.")

    if n2 is not None:
        n2 = int(n2)
        if n2 <= 0:
            raise ValueError("n2 must be a positive integer when provided.")

    a1_lo, a1_hi = _parse_angle_range(a1, "a1", np.pi)
    a2_lo, a2_hi = _parse_angle_range(a2, "a2", 2.0 * np.pi)

    a1_range = a1_hi - a1_lo
    a2_range = a2_hi - a2_lo

    da1 = a1_range / n1
    da2 = da1 if (n2 is None) else (a2_range / n2)

    # -----------------------------
    # Reproduce MATLAB loop over i:
    # for i = (-pi/2 + a1(1)) : da1 : (a1(2) - pi/2)
    # Using linspace gives n1+1 samples including endpoints.
    # -----------------------------
    i_vals = np.linspace(-np.pi / 2.0 + a1_lo, a1_hi - np.pi / 2.0, n1 + 1)

    x_list = []
    y_list = []
    z_list = []

    for i in i_vals:
        # n = ceil(a2ran*cos(i)/da2)
        # (MATLAB: n=ceil(a2ran*cos(i)/da2);)
        n = int(np.ceil(a2_range * np.cos(i) / da2)) if da2 != 0 else 0

        if n == 0:
            # MATLAB special case:
            # x=[x 0]; y=[y 0]; z=[z 1];
            x_list.append(0.0)
            y_list.append(0.0)
            z_list.append(1.0)
        else:
            # MATLAB:
            # phi = a2(1):a2ran/n:a2(2)  -> n+1 samples incl endpoints
            phi = np.linspace(a2_lo, a2_hi, n + 1)

            cos_i = np.cos(i)
            sin_i = np.sin(i)

            # x = sin(phi)*cos(i)
            x_ring = np.sin(phi) * cos_i

            # y = cos(phi)*(-cos(i))
            y_ring = np.cos(phi) * (-cos_i)

            # z = sin(i) repeated
            z_ring = np.full(phi.shape, sin_i, dtype=float)

            x_list.extend(x_ring.tolist())
            y_list.extend(y_ring.tolist())
            z_list.extend(z_ring.tolist())

    xyz = np.column_stack((np.asarray(x_list, dtype=float),
                           np.asarray(y_list, dtype=float),
                           np.asarray(z_list, dtype=float)))
    return xyz


ArrayLike = Union[np.ndarray, list, tuple]

def gentang(
    v: ArrayLike,
    angle_span: float = 2 * np.pi,
    n_systems: int = 1,
    start_angle: float = 0.0,
    *,
    squeeze: bool = True,
    return_matrices: bool = False,
    eps: float = 1e-12,
) -> np.ndarray:
    """
    Generate (for each input vector) an orthonormal coordinate system based on
    a circle of tangent directions in the plane perpendicular to the vector.

    Parameters
    ----------
    v : array-like, shape (n, 3) or (3,)
        Input vectors. They do NOT need to be unit length (will be normalized).
    angle_span : float, default 2*pi
        Total angular span (radians) over which to distribute tangent directions.
        Like MATLAB: points are uniform on [start_angle, start_angle + angle_span).
    n_systems : int, default 1
        Number of coordinate systems to generate per input vector.
        (Equivalent to MATLAB 'an'.) With n_systems=1 you get exactly one basis per v.
    start_angle : float, default 0
        Starting angle (radians) for the tangent sweep.
    squeeze : bool, default True
        If True and n_systems==1, return shape (n, 9) (or (n, 3, 3) if return_matrices).
        Otherwise return shape (n, n_systems, 9) (or (n, n_systems, 3, 3)).
    return_matrices : bool, default False
        If False: returns direction cosines flattened as [Vn, Vs1, Vs2] (length 9).
        If True: returns 3x3 basis matrices whose columns are [Vs1, Vs2, Vn].
    eps : float, default 1e-12
        Threshold for handling near-singular cases (vectors nearly aligned with X axis).

    Returns
    -------
    out : ndarray
        If return_matrices=False:
            - (n, 9) if squeeze and n_systems==1
            - (n, n_systems, 9) otherwise
          with columns [lN mN nN lS1 mS1 nS1 lS2 mS2 nS2]
        If return_matrices=True:
            - (n, 3, 3) if squeeze and n_systems==1
            - (n, n_systems, 3, 3) otherwise
          where columns are [Vn, Vs1, Vs2].
    """
    v = np.asarray(v, dtype=float)
    if v.ndim == 1:
        if v.shape[0] != 3:
            raise ValueError("If v is 1D, it must have length 3.")
        v = v.reshape(1, 3)
    if v.ndim != 2 or v.shape[1] != 3:
        raise ValueError("v must have shape (n, 3) or (3,).")

    if n_systems < 1:
        raise ValueError("n_systems must be >= 1.")

    # Normalize input vectors (Vn)
    norms = np.linalg.norm(v, axis=1)
    if np.any(norms < eps):
        raise ValueError("v contains (near-)zero vectors; cannot normalize.")
    vn = v / norms[:, None]  # (n,3)

    l = vn[:, 0]
    m = vn[:, 1]
    n_ = vn[:, 2]

    # Angles like MATLAB: uniform points on [start_angle, start_angle + angle_span)
    angles = start_angle + (np.arange(n_systems) * (angle_span / n_systems))
    y = np.cos(angles)  # (n_systems,)
    z = np.sin(angles)  # (n_systems,)

    # Broadcast to (n, n_systems)
    Y = y[None, :]
    Z = z[None, :]

    s = np.sqrt(m**2 + n_**2)  # (n,)
    S = s[:, None]             # (n,1)

    # Tangent direction Vs1 (matches MATLAB's csize=3 "xp,yp,zp" construction)
    # Special case when vn is aligned with X axis (m^2+n^2 ~ 0)
    vs1 = np.empty((vn.shape[0], n_systems, 3), dtype=float)

    rows_singular = s < eps
    rows_regular = ~rows_singular

    if np.any(rows_regular):
        lr = l[rows_regular][:, None]
        mr = m[rows_regular][:, None]
        nr = n_[rows_regular][:, None]
        Sr = S[rows_regular]

        xp = -Z * Sr
        yp = (-Y * nr + Z * lr * mr) / Sr
        zp = ( Y * mr + Z * lr * nr) / Sr
        vs1[rows_regular, :, 0] = xp
        vs1[rows_regular, :, 1] = yp
        vs1[rows_regular, :, 2] = zp

    if np.any(rows_singular):
        # vn is (approximately) ±X; pick a clean circle in YZ plane
        vs1[rows_singular, :, 0] = 0.0
        vs1[rows_singular, :, 1] = y
        vs1[rows_singular, :, 2] = z

    # Vs2 from cross product to guarantee orthonormal, right-handed basis:
    # columns [Vn, Vs1, Vs2] with Vs2 = Vn x Vs1 => Vs1 x Vs2 = Vn
    vn_rep = np.repeat(vn[:, None, :], n_systems, axis=1)  # (n, n_systems, 3)
    vs2 = np.cross(vn_rep, vs1)

    # Normalize Vs2 defensively (numerical safety)
    vs2_norm = np.linalg.norm(vs2, axis=2, keepdims=True)
    if np.any(vs2_norm < eps):
        raise ValueError("Failed to construct Vs2 for some vectors (degenerate cross product).")
    vs2 = vs2 / vs2_norm

    if return_matrices:
        # Basis matrices with columns [X, Y, Z] = [Vs1, Vs2, Vn]
        # This is right-handed because Vs2 = Vn x Vs1  =>  Vs1 x Vs2 = Vn
        basis = np.stack([vs1, vs2, vn_rep], axis=-1)  # (n, n_systems, 3, 3)
        if squeeze and n_systems == 1:
            return basis[:, 0, :, :]
        return basis

    # Flatten like MATLAB: [lN mN nN lS1 mS1 nS1 lS2 mS2 nS2]
    out = np.concatenate([vn_rep, vs1, vs2], axis=2)  # (n, n_systems, 9)

    if squeeze and n_systems == 1:
        return out[:, 0, :]
    return out

def build_R_list(*, n1: int = 20, n2: int = 40) -> np.ndarray:
    """
    Returns an array of shape (N, 3, 3) of orthonormal rotation matrices.
    Convention: the 3rd column of R is the plane normal (local z-axis).

    This is robust to different gentang() return conventions.
    """
    xyz = dcmesh_xyz(a1=np.pi, n1=n1, a2=2*np.pi, n2=n2)
    xyz = np.unique(np.round(xyz, 12), axis=0)

    g = gentang(xyz, return_matrices=True)

    # Case A: gentang returns R_list directly (like your MAIN.py usage)
    if isinstance(g, np.ndarray):
        R_list = g

    # Case B: gentang returns a tuple/list of outputs
    elif isinstance(g, (tuple, list)):
        # Common patterns:
        #   (R_list,)
        #   (tang, R_list)
        #   (xyz2, tang, R_list)
        # We take the last item that looks like (N,3,3).
        R_list = None
        for item in reversed(g):
            arr = np.asarray(item)
            if arr.ndim == 3 and arr.shape[1:] == (3, 3):
                R_list = arr
                break
        if R_list is None:
            raise ValueError(
                "gentang(return_matrices=True) did not return an array shaped (N,3,3). "
                f"Got types/shapes: {[getattr(np.asarray(x), 'shape', None) for x in g]}"
            )
    else:
        raise TypeError(f"Unexpected gentang return type: {type(g)}")

    return np.asarray(R_list, dtype=float)
