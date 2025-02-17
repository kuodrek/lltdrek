import numpy as np
import math
from typing import List, Union, Optional
from dataclasses import dataclass, field
from src.models.flight_condition import FlightCondition
import src.utils.geometry as geo


@dataclass(repr=False, eq=False, match_args=False, slots=True)
class Wing:
    spans: List[float]
    chords: List[float]
    offsets: List[float]
    twist_angles: List[float]
    dihedral_angles: List[float]
    airfoils: List[str]
    surface_name: str
    N_panels: Union[int, List[int]]
    x_pos: float = 0
    z_pos: float = 0
    distribution_type: str = 'linear'
    sweep_check: bool = False
    # Properties obtained through get_mesh method
    total_span: float = field(init=False)
    span_panel_numbers: List[int] = field(init=False)
    collocation_points: np.ndarray = field(init=False)
    vertice_points: np.ndarray = field(init=False)
    u_a: np.ndarray = field(init=False)
    u_n: np.ndarray = field(init=False)
    u_s: np.ndarray = field(init=False)
    cp_lengths: np.ndarray = field(init=False)
    cp_dsl: np.ndarray = field(init=False)
    cp_areas: np.ndarray = field(init=False)
    cp_chords: np.ndarray = field(init=False)
    cp_macs: np.ndarray = field(init=False)
    cp_reynolds: np.ndarray = field(init=False)
    cp_airfoils: List = field(init=False)
    airfoil_data: dict = field(init=False)
    MAC: float = field(init=False)
    partition_areas: np.ndarray = field(init=False)
    total_area: float = field(init=False)
    AR: float = field(init=False)
    parent_wing: Optional[str] = field(init=False)

    allowed_distribution_types = [
        "linear",
        "cosine"
    ]

    def __post_init__(self) -> None:
        if self.distribution_type not in self.allowed_distribution_types:
            raise ValueError("Invalid 'distribution_type' value.")
        # Convert degree to rad
        self.twist_angles = [angle * np.pi / 180 for angle in self.twist_angles]
        self.dihedral_angles = [angle * np.pi / 180 for angle in self.dihedral_angles]
        self.parent_wing = None # Usado para asas espelhadas
    

    def __repr__(self) -> None:
        return self.surface_name


    def generate_mesh(self) -> None:
        self.total_span = sum(self.spans)
        # Number of panels of each partition
        self.span_panel_numbers = [math.ceil(value / self.total_span * self.N_panels) for value in self.spans]
        # A quantidade de paineis da última partição é definida como N_panels - (paineis de todas as outras partições)
        self.span_panel_numbers[-1] = self.N_panels - sum(self.span_panel_numbers) + self.span_panel_numbers[-1]

        self.partition_areas = np.zeros(len(self.spans))
        
        self.collocation_points = np.zeros([self.N_panels, 3])
        self.vertice_points = np.zeros([self.N_panels + 1, 3])

        self.u_a = np.zeros([self.N_panels, 3]) # Vetor unitário colinear à corda
        self.u_n = np.zeros([self.N_panels, 3]) # Vetor unitário normal à corda
        self.u_s = np.zeros([self.N_panels, 3]) # Vetor perpendicular ao plano do perfil

        self.cp_lengths = np.zeros([self.N_panels, 3]) # Vetor de comprimento de cada painel
        self.cp_dsl = np.zeros([self.N_panels, 3]) # Vetor de comprimento ao longo da envergadura adimensional (dimensionless spanwise length vector)

        self.cp_areas = np.zeros(self.N_panels) # Área de cada painel
        self.cp_chords = np.zeros(self.N_panels) # Distribuição de cordas dos pontos de colocação
        self.cp_macs = np.zeros(self.N_panels) # Corda média aerodinâmica de cada painel

        # TODO: explain what these variables mean
        idx_n = 0
        span_incremental = 0
        height_incremental = 0
        MAC = 0
        for i, span_partition in enumerate(self.spans):
            n = self.span_panel_numbers[i]
            chord_i = self.chords[i]
            chord_ii = self.chords[i+1]
            offset_i = self.offsets[i]
            offset_ii = self.offsets[i+1]
            twist_i = self.twist_angles[i]
            twist_ii = self.twist_angles[i+1]
            dihedral_partition = self.dihedral_angles[i]

            # Ângulo de enflechamento [rad] da partição atual
            sweep_partition = np.arctan( (0.25*chord_ii - 0.25*chord_i + offset_i)/span_partition ) if self.sweep_check is True else 0

            span_y_distr = geo.get_span_y_distr(n, span_partition, self.distribution_type)
            cp_y_component = span_y_distr['collocation_points']
            vp_y_component = span_y_distr['vertice_points']

            span_x_distr = geo.get_span_x_distr(cp_y_component,vp_y_component,chord_i,chord_ii,offset_i,offset_ii,span_partition)
            cp_x_component = span_x_distr['collocation_points']
            vp_x_component = span_x_distr['vertice_points']

            span_z_distr = geo.get_span_z_distr(cp_y_component,vp_y_component,dihedral_partition)
            cp_z_component = span_z_distr['collocation_points']
            vp_z_component = span_z_distr['vertice_points']

            cp_twist_distr = geo.get_twist_distr(cp_y_component, twist_i, twist_ii, span_partition)

            # Alocar o primeiro ponto em vertice_points
            if i == 0:
                self.vertice_points[0][0] = vp_x_component[0]
                self.vertice_points[0][1] = vp_y_component[0]
                self.vertice_points[0][2] = vp_z_component[0]

            for j, _ in enumerate(cp_y_component):
                # Componente Y dos CP / VP
                self.vertice_points[idx_n+j+1][1] = span_incremental + vp_y_component[j+1]
                self.collocation_points[idx_n+j][1] = span_incremental + cp_y_component[j]
                # Componente X
                self.vertice_points[idx_n+j+1][0] = vp_x_component[j+1]
                self.collocation_points[idx_n+j][0] = cp_x_component[j]
                # Componente Z
                self.vertice_points[idx_n+j+1][2] = height_incremental + vp_z_component[j+1]
                self.collocation_points[idx_n+j][2] = height_incremental + cp_z_component[j]
                # Ângulo de torção (twist geométrico) do ponto de colocação [testar]
                twist_cp = cp_twist_distr[j]

                # Vetores de posicionamento do painel (colinear e perpendicular à corda)
                euler_matrix = geo.get_euler_matrix(dihedral_partition, twist_cp, sweep_partition)
                self.u_a[idx_n+j] = euler_matrix.dot(np.array([1, 0, 0]))
                self.u_n[idx_n+j] = euler_matrix.dot(np.array([0, 0, 1]))
                self.u_s[idx_n+j] = np.cross(self.u_a[idx_n+j], self.u_n[idx_n+j])

                chord_vp_j = geo.get_local_chord(vp_y_component[j], chord_i, chord_ii, span_partition)
                chord_vp_jj = geo.get_local_chord(vp_y_component[j+1], chord_i, chord_ii, span_partition)
                chord_cp = geo.get_local_chord(cp_y_component[j], chord_i, chord_ii, span_partition)
                self.cp_chords[idx_n+j] = chord_cp

                # Corda média aerodinâmica  do painel
                cp_mac = (2/3)*(chord_vp_j ** 2 + chord_vp_j * chord_vp_jj + chord_vp_jj ** 2)/(chord_vp_j + chord_vp_jj)
                self.cp_macs[idx_n+j] = cp_mac
                for k in range(3):
                    self.cp_lengths[idx_n+j][k] = self.vertice_points[idx_n+j+1][k] - self.vertice_points[idx_n+j][k]

                # Vetor comprimento do painel
                cp_length_x = self.cp_lengths[idx_n+j][0]
                cp_length_y = self.cp_lengths[idx_n+j][1]
                cp_length_z = self.cp_lengths[idx_n+j][2]

                # Área do painel
                cp_area = 0.5*(chord_vp_j + chord_vp_jj)*math.sqrt( cp_length_y ** 2 + cp_length_z ** 2 )
                self.cp_areas[idx_n+j] = cp_area
                self.partition_areas[i] += cp_area

                # Cálculo do dimensionless spanwise length
                for k in range(3): 
                    self.cp_dsl[idx_n+j][k] = (cp_mac * self.cp_lengths[idx_n+j][k]) / cp_area

            # Razão de afilamento da partição atual
            partition_lambda = chord_ii / chord_i
            # Corda média aerodinâmica da asa
            MAC += self.partition_areas[i]*(2/3)*chord_i*(1+partition_lambda+partition_lambda**2)/(1+partition_lambda)
            idx_n += n
            span_incremental += span_partition
            height_incremental += vp_z_component[-1]
    
        # Add coordinates into cp's and vp's
        self.collocation_points[:,0] += self.x_pos
        self.vertice_points[:,0] += self.x_pos
        self.collocation_points[:,2] += self.z_pos
        self.vertice_points[:,2] += self.z_pos
        
        # Global values
        self.total_area = sum(self.partition_areas)
        self.MAC = MAC / sum(self.partition_areas)
        self.AR = (2 * self.total_span) ** 2 / (2 * sum(self.partition_areas))


    def setup_airfoil_data(self, flight_condition: FlightCondition, airfoils_data_dict: dict) -> None:
        """
        TODO: Explicar este método e suas variáveis
        """
        self.cp_reynolds = np.zeros(self.N_panels)
        self.cp_airfoils = []
        idx_n = 0
        for i, span_partition in enumerate(self.spans):
            n = self.span_panel_numbers[i]
            airfoil_i = self.airfoils[i]
            airfoil_ii = self.airfoils[i+1]
            for j in range(self.span_panel_numbers[i]):
                if airfoil_i == airfoil_ii:
                    self.cp_airfoils.append([1, [airfoil_i]])
                else:
                    merge_parameter = (self.collocation_points[idx_n+j][1]-self.vertice_points[idx_n][1]) / span_partition
                    self.cp_airfoils.append([merge_parameter, [airfoil_i, airfoil_ii]])
            idx_n += n

        for i, panel_chord in enumerate(self.cp_chords):
            reynolds_number = panel_chord * flight_condition.V_inf / flight_condition.nu
            self.cp_reynolds[i] = reynolds_number
        
        unique_airfoils = list(set(self.airfoils))
        self.airfoil_data = {}
        for airfoil, airfoil_data in airfoils_data_dict.items():
            if airfoil in unique_airfoils:
                self.airfoil_data[airfoil] = airfoil_data
