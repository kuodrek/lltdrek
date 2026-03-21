from abc import ABC, abstractmethod

import numpy as np
import numpy.linalg as npla

# ---------------------------------------------------------------------------
# Private helpers — shared by all concrete implementations
# ---------------------------------------------------------------------------


def _biot_savart(
    cp: np.ndarray,
    vp1: np.ndarray,
    vp2: np.ndarray,
    mac: float,
    v_inf: np.ndarray,
) -> tuple[np.ndarray, bool]:
    """Pure Biot-Savart formula for one horseshoe vortex panel.

    Computes the velocity induced at collocation point `cp` by the horseshoe
    vortex whose bound segment runs from `vp1` to `vp2`, with two semi-infinite
    trailing legs aligned with `v_inf`.

    Returns
    -------
    velocity : (3,) ndarray
    same_panel_check : bool
        True when the bound-vortex denominator is ~0 (self-induction panel).
        The bound vortex contribution is omitted in that case.
    """
    ri1j = cp - vp1
    ri2j = cp - vp2

    ri1j_abs = npla.norm(ri1j)
    ri2j_abs = npla.norm(ri2j)

    r1_cross = np.cross(v_inf, ri1j)
    r2_cross = np.cross(v_inf, ri2j)
    r1_dot = np.dot(v_inf, ri1j)
    r2_dot = np.dot(v_inf, ri2j)

    velocity = (
        mac / (4 * np.pi) * (r2_cross / (ri2j_abs * (ri2j_abs - r2_dot)) - r1_cross / (ri1j_abs * (ri1j_abs - r1_dot)))
    )

    r12_cross = np.cross(ri1j, ri2j)
    r12_dot = np.dot(ri1j, ri2j)
    bound_den = ri1j_abs * ri2j_abs * (ri1j_abs * ri2j_abs + r12_dot)

    if not np.isclose(bound_den, 0, atol=1e-10):
        velocity += mac / (4 * np.pi) * (ri1j_abs + ri2j_abs) * r12_cross / bound_den
        same_panel_check = False
    else:
        same_panel_check = True

    return velocity, same_panel_check


def _reflect_ground(
    vp1: np.ndarray,
    vp2: np.ndarray,
    v_inf: np.ndarray,
    h: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reflect a vortex panel below the ground plane (image method).

    Transforms vertex coordinates and freestream direction for the image vortex:
    - Image z-coordinate: z_img = 2h + 2*z_real  (matches legacy formula)
    - v_inf z-component is negated

    Returns copies — does not modify inputs.
    """
    vp1_img = vp1.copy()
    vp2_img = vp2.copy()
    vp1_img[2] = 2 * h + 2 * vp1[2]
    vp2_img[2] = 2 * h + 2 * vp2[2]

    v_inf_img = v_inf.copy()
    v_inf_img[2] = -v_inf_img[2]

    return vp1_img, vp2_img, v_inf_img


# ---------------------------------------------------------------------------
# ABC
# ---------------------------------------------------------------------------


class VelocityCalculator(ABC):
    @abstractmethod
    def get_induced_velocity_distribution(
        self,
        collocation_points: np.ndarray,
        cp_macs: np.ndarray,
        vertice_points: np.ndarray,
        freestream_velocity: np.ndarray,
        is_mirrored: bool = False,
        ground_effect: bool = False,
        h: float = 0.0,
    ) -> np.ndarray:
        """Induced velocity matrix from all vortex panels on all collocation points.

        Parameters
        ----------
        collocation_points : (N_cp, 3)
        cp_macs : (N_cp,)
        vertice_points : (N_cp+1, 3)
        freestream_velocity : (N_cp, 3)  — one vector per collocation point
        is_mirrored : bool
            True when the surface is a mirrored copy; reverses panel ordering so
            vortex direction is consistent with the original surface.
        ground_effect : bool
            Adds image-vortex contribution via the ground reflection method.
        h : float
            Height parameter for ground effect (unused when ground_effect=False).

        Returns
        -------
        v_ij : (N_cp, N_cp, 3)
            v_ij[i, j] is the velocity induced at cp_i by panel j.
        """


# ---------------------------------------------------------------------------
# Concrete implementations
# ---------------------------------------------------------------------------


class LoopsVelocityCalculator(VelocityCalculator):
    """Induced velocity distribution computed with nested Python for-loops.

    Direct port of the legacy ``get_induced_velocity_distribution``.
    Validated baseline — prefer this when correctness matters over speed.
    """

    def get_induced_velocity_distribution(
        self,
        collocation_points: np.ndarray,
        cp_macs: np.ndarray,
        vertice_points: np.ndarray,
        freestream_velocity: np.ndarray,
        is_mirrored: bool = False,
        ground_effect: bool = False,
        h: float = 0.0,
    ) -> np.ndarray:
        n_cp = len(collocation_points)
        n_panels = len(vertice_points) - 1
        result = np.zeros((n_cp, n_panels, 3))

        for i in range(n_cp):
            cp_i = collocation_points[i]
            mac_i = cp_macs[i]
            v_inf_i = freestream_velocity[i]

            for j in range(n_panels):
                if is_mirrored:
                    vp_j = vertice_points[j + 1]
                    vp_jj = vertice_points[j]
                else:
                    vp_j = vertice_points[j]
                    vp_jj = vertice_points[j + 1]

                vel, _ = _biot_savart(cp_i, vp_j, vp_jj, mac_i, v_inf_i)

                if ground_effect:
                    vp_j_img, vp_jj_img, v_inf_img = _reflect_ground(vp_j, vp_jj, v_inf_i, h)
                    vel += _biot_savart(cp_i, vp_j_img, vp_jj_img, mac_i, v_inf_img)[0]

                result[i, j] = vel

        return result


class NumpyVelocityCalculator(VelocityCalculator):
    """Induced velocity distribution using fully vectorized numpy operations.

    Implements the same equations as ``LoopsVelocityCalculator`` with no
    Python loops over panels. Validate against ``LoopsVelocityCalculator``
    before use.
    """

    def get_induced_velocity_distribution(
        self,
        collocation_points: np.ndarray,
        cp_macs: np.ndarray,
        vertice_points: np.ndarray,
        freestream_velocity: np.ndarray,
        is_mirrored: bool = False,
        ground_effect: bool = False,
        h: float = 0.0,
    ) -> np.ndarray:
        raise NotImplementedError
