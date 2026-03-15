import math
from typing import List, Optional, Union

import numpy as np

import lltdrek.utils.geometry as geo
from new.geometry.wing import Wing


class WingLLT(Wing):
    """Lifting Line Theory wing implementation.

    Builds a panel mesh using LLT conventions: horseshoe vortex panels
    distributed spanwise with cosine or linear spacing. Euler angle matrices
    (dihedral, twist, sweep) define each panel's orientation vectors.
    """

    _ALLOWED_DISTRIBUTION_TYPES = ("linear", "cosine")

    def __init__(
        self,
        spans: List[float],
        chords: List[float],
        offsets: List[float],
        twist_angles: List[float],
        dihedral_angles: List[float],
        airfoils: List[str],
        surface_name: str,
        N_panels: Union[int, List[int]],
        x_pos: float = 0,
        z_pos: float = 0,
        distribution_type: str = "linear",
        sweep_check: bool = False,
    ) -> None:
        if distribution_type not in self._ALLOWED_DISTRIBUTION_TYPES:
            raise ValueError(f"Invalid 'distribution_type'. Choose from {self._ALLOWED_DISTRIBUTION_TYPES}.")

        self.spans = spans
        self.chords = chords
        self.offsets = offsets
        self.airfoils = airfoils
        self._surface_name = surface_name
        self._N_panels = N_panels
        self._x_pos = x_pos
        self.z_pos = z_pos
        self.distribution_type = distribution_type
        self.sweep_check = sweep_check

        # Convert degrees to radians
        self.twist_angles = [a * np.pi / 180 for a in twist_angles]
        self.dihedral_angles = [a * np.pi / 180 for a in dihedral_angles]

        # Set after generate_mesh()
        self._total_span: float = None
        self._total_area: float = None
        self._MAC: float = None
        self._AR: float = None
        self._collocation_points: np.ndarray = None
        self._vertice_points: np.ndarray = None
        self._u_a: np.ndarray = None
        self._u_n: np.ndarray = None
        self._u_s: np.ndarray = None
        self._cp_lengths: np.ndarray = None
        self._cp_dsl: np.ndarray = None
        self._cp_areas: np.ndarray = None
        self._cp_chords: np.ndarray = None
        self._cp_macs: np.ndarray = None
        self._partition_areas: np.ndarray = None
        self._span_panel_numbers: List[int] = None

        # Set after setup_airfoil_data()
        self._cp_reynolds: np.ndarray = None
        self._cp_airfoils: list = None

        # Set by WingPool when creating mirrored copies
        self._parent_wing: Optional[str] = None

        self.generate_mesh()

    def __repr__(self) -> str:
        return self._surface_name

    # -------------------------------------------------------------------------
    # Wing metadata properties
    # -------------------------------------------------------------------------

    @property
    def surface_name(self) -> str:
        return self._surface_name

    @property
    def N_panels(self) -> int:
        return self._N_panels

    @property
    def x_pos(self) -> float:
        return self._x_pos

    # -------------------------------------------------------------------------
    # Computed wing-level metrics
    # -------------------------------------------------------------------------

    @property
    def total_span(self) -> float:
        return self._total_span

    @property
    def total_area(self) -> float:
        return self._total_area

    @property
    def MAC(self) -> float:
        return self._MAC

    @property
    def AR(self) -> float:
        return self._AR

    # -------------------------------------------------------------------------
    # Panel mesh outputs (read-write — WingPool mutates mirrored copies)
    # -------------------------------------------------------------------------

    @property
    def collocation_points(self) -> np.ndarray:
        return self._collocation_points

    @collocation_points.setter
    def collocation_points(self, value: np.ndarray) -> None:
        self._collocation_points = value

    @property
    def vertice_points(self) -> np.ndarray:
        return self._vertice_points

    @vertice_points.setter
    def vertice_points(self, value: np.ndarray) -> None:
        self._vertice_points = value

    @property
    def u_a(self) -> np.ndarray:
        return self._u_a

    @u_a.setter
    def u_a(self, value: np.ndarray) -> None:
        self._u_a = value

    @property
    def u_n(self) -> np.ndarray:
        return self._u_n

    @u_n.setter
    def u_n(self, value: np.ndarray) -> None:
        self._u_n = value

    @property
    def u_s(self) -> np.ndarray:
        return self._u_s

    @u_s.setter
    def u_s(self, value: np.ndarray) -> None:
        self._u_s = value

    @property
    def cp_lengths(self) -> np.ndarray:
        return self._cp_lengths

    @cp_lengths.setter
    def cp_lengths(self, value: np.ndarray) -> None:
        self._cp_lengths = value

    @property
    def cp_dsl(self) -> np.ndarray:
        return self._cp_dsl

    @cp_dsl.setter
    def cp_dsl(self, value: np.ndarray) -> None:
        self._cp_dsl = value

    @property
    def cp_areas(self) -> np.ndarray:
        return self._cp_areas

    @property
    def cp_chords(self) -> np.ndarray:
        return self._cp_chords

    @property
    def cp_macs(self) -> np.ndarray:
        return self._cp_macs

    # -------------------------------------------------------------------------
    # Airfoil data properties
    # -------------------------------------------------------------------------

    @property
    def cp_reynolds(self) -> np.ndarray:
        return self._cp_reynolds

    @property
    def cp_airfoils(self) -> list:
        return self._cp_airfoils

    # -------------------------------------------------------------------------
    # Cross-wing references
    # -------------------------------------------------------------------------

    @property
    def parent_wing(self) -> Optional[str]:
        return self._parent_wing

    @parent_wing.setter
    def parent_wing(self, value: Optional[str]) -> None:
        self._parent_wing = value

    # -------------------------------------------------------------------------
    # Lifecycle methods
    # -------------------------------------------------------------------------

    def generate_mesh(self) -> None:
        """Build the LLT panel mesh for this wing surface."""
        self._total_span = sum(self.spans)

        # Number of panels per partition (proportional to span fraction)
        self._span_panel_numbers = [math.ceil(s / self._total_span * self._N_panels) for s in self.spans]
        # Last partition absorbs rounding remainder
        self._span_panel_numbers[-1] = self._N_panels - sum(self._span_panel_numbers) + self._span_panel_numbers[-1]

        self._partition_areas = np.zeros(len(self.spans))
        self._collocation_points = np.zeros([self._N_panels, 3])
        self._vertice_points = np.zeros([self._N_panels + 1, 3])

        self._u_a = np.zeros([self._N_panels, 3])
        self._u_n = np.zeros([self._N_panels, 3])
        self._u_s = np.zeros([self._N_panels, 3])

        self._cp_lengths = np.zeros([self._N_panels, 3])
        self._cp_dsl = np.zeros([self._N_panels, 3])
        self._cp_areas = np.zeros(self._N_panels)
        self._cp_chords = np.zeros(self._N_panels)
        self._cp_macs = np.zeros(self._N_panels)

        idx_n = 0
        span_incremental = 0
        height_incremental = 0
        MAC = 0

        for i, span_partition in enumerate(self.spans):
            n = self._span_panel_numbers[i]
            chord_i = self.chords[i]
            chord_ii = self.chords[i + 1]
            offset_i = self.offsets[i]
            offset_ii = self.offsets[i + 1]
            twist_i = self.twist_angles[i]
            twist_ii = self.twist_angles[i + 1]
            dihedral_partition = self.dihedral_angles[i]

            sweep_partition = (
                np.arctan((0.25 * chord_ii - 0.25 * chord_i + offset_i) / span_partition) if self.sweep_check else 0
            )

            span_y_distr = geo.get_span_y_distr(n, span_partition, self.distribution_type)
            cp_y = span_y_distr["collocation_points"]
            vp_y = span_y_distr["vertice_points"]

            span_x_distr = geo.get_span_x_distr(cp_y, vp_y, chord_i, chord_ii, offset_i, offset_ii, span_partition)
            cp_x = span_x_distr["collocation_points"]
            vp_x = span_x_distr["vertice_points"]

            span_z_distr = geo.get_span_z_distr(cp_y, vp_y, dihedral_partition)
            cp_z = span_z_distr["collocation_points"]
            vp_z = span_z_distr["vertice_points"]

            cp_twist_distr = geo.get_twist_distr(cp_y, twist_i, twist_ii, span_partition)

            if i == 0:
                self._vertice_points[0] = [vp_x[0], vp_y[0], vp_z[0]]

            for j in range(len(cp_y)):
                self._vertice_points[idx_n + j + 1] = [
                    vp_x[j + 1],
                    span_incremental + vp_y[j + 1],
                    height_incremental + vp_z[j + 1],
                ]
                self._collocation_points[idx_n + j] = [
                    cp_x[j],
                    span_incremental + cp_y[j],
                    height_incremental + cp_z[j],
                ]

                euler_matrix = geo.get_euler_matrix(dihedral_partition, cp_twist_distr[j], sweep_partition)
                self._u_a[idx_n + j] = euler_matrix.dot(np.array([1, 0, 0]))
                self._u_n[idx_n + j] = euler_matrix.dot(np.array([0, 0, 1]))
                self._u_s[idx_n + j] = np.cross(self._u_a[idx_n + j], self._u_n[idx_n + j])

                chord_vp_j = geo.get_local_chord(vp_y[j], chord_i, chord_ii, span_partition)
                chord_vp_jj = geo.get_local_chord(vp_y[j + 1], chord_i, chord_ii, span_partition)
                chord_cp = geo.get_local_chord(cp_y[j], chord_i, chord_ii, span_partition)
                self._cp_chords[idx_n + j] = chord_cp

                cp_mac = (
                    (2 / 3) * (chord_vp_j**2 + chord_vp_j * chord_vp_jj + chord_vp_jj**2) / (chord_vp_j + chord_vp_jj)
                )
                self._cp_macs[idx_n + j] = cp_mac

                self._cp_lengths[idx_n + j] = self._vertice_points[idx_n + j + 1] - self._vertice_points[idx_n + j]

                cp_length_y = self._cp_lengths[idx_n + j][1]
                cp_length_z = self._cp_lengths[idx_n + j][2]
                cp_area = 0.5 * (chord_vp_j + chord_vp_jj) * math.sqrt(cp_length_y**2 + cp_length_z**2)
                self._cp_areas[idx_n + j] = cp_area
                self._partition_areas[i] += cp_area

                self._cp_dsl[idx_n + j] = (cp_mac * self._cp_lengths[idx_n + j]) / cp_area

            partition_lambda = chord_ii / chord_i
            MAC += (
                self._partition_areas[i]
                * (2 / 3)
                * chord_i
                * (1 + partition_lambda + partition_lambda**2)
                / (1 + partition_lambda)
            )
            idx_n += n
            span_incremental += span_partition
            height_incremental += vp_z[-1]

        self._collocation_points[:, 0] += self._x_pos
        self._vertice_points[:, 0] += self._x_pos
        self._collocation_points[:, 2] += self.z_pos
        self._vertice_points[:, 2] += self.z_pos

        self._total_area = float(sum(self._partition_areas))
        self._MAC = MAC / self._total_area
        self._AR = (2 * self._total_span) ** 2 / (2 * self._total_area)

    def _apply_flight_condition(self, flight_condition) -> None:
        """Assign per-panel airfoils and compute Reynolds numbers from the flight condition."""
        self._cp_reynolds = np.zeros(self._N_panels)
        self._cp_airfoils = []

        idx_n = 0
        for i, span_partition in enumerate(self.spans):
            airfoil_i = self.airfoils[i]
            airfoil_ii = self.airfoils[i + 1]
            for j in range(self._span_panel_numbers[i]):
                if airfoil_i == airfoil_ii:
                    self._cp_airfoils.append([1, [airfoil_i]])
                else:
                    merge_parameter = (
                        self._collocation_points[idx_n + j][1] - self._vertice_points[idx_n][1]
                    ) / span_partition
                    self._cp_airfoils.append([merge_parameter, [airfoil_i, airfoil_ii]])
            idx_n += self._span_panel_numbers[i]

        for i, panel_chord in enumerate(self._cp_chords):
            self._cp_reynolds[i] = panel_chord * flight_condition.V_inf / flight_condition.nu
