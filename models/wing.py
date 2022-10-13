import numpy as np
from models.flight_condition import FlightCondition
import utils.geometry as geo
import math
from typing import List, Union
from dataclasses import dataclass, field


@dataclass
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
    distribution_type: str = 'Linear'
    sweep_check: bool = False
    # Properties obtained through get_mesh method
    total_span: float = field(init=False)
    span_panel_numbers: List[int] = field(init=False)
    collocation_points: np.ndarray = field(init=False)
    vertice_points: np.ndarray = field(init=False)
    u_a: np.ndarray = field(init=False)
    u_n: np.ndarray = field(init=False)
    cp_lengths: np.ndarray = field(init=False)
    cp_dsl: np.ndarray = field(init=False)
    cp_areas: np.ndarray = field(init=False)
    cp_chords: np.ndarray = field(init=False)
    cp_macs: np.ndarray = field(init=False)
    cp_reynolds: np.ndarray = field(init=False)
    cp_airfoils: np.ndarray = field(init=False)
    MAC: float = field(init=False)
    partition_areas: np.ndarray = field(init=False)
    total_area: float = field(init=False)
    AR: float = field(init=False)

    def __post_init__(self):
        # Convert degree to rad
        self.twist_angles = [angle * np.pi / 180 for angle in self.twist_angles]
        self.dihedral_angles = [angle * np.pi / 180 for angle in self.dihedral_angles]
    
    def generate_mesh(self):
        N_panels = self.N_panels
        spans = self.spans
        chords = self.chords
        offsets = self.offsets
        twist_angles = self.twist_angles
        dihedral_angles = self.dihedral_angles
        airfoils = self.airfoils
        distribution_type = self.distribution_type
        sweep_check = self.sweep_check

        total_span = sum(spans)
        # Number of panels of each partition
        span_panel_numbers = [value / total_span * N_panels for value in spans]
        span_panel_numbers = [math.ceil(int(i)) for i in span_panel_numbers]
        # Add a panel to compensate rounding
        if sum(span_panel_numbers) == N_panels - 1:
            span_panel_numbers[-1] += 1

        partition_areas = np.zeros(len(spans))
        
        collocation_points = np.zeros([N_panels, 3])
        vertice_points = np.zeros([N_panels + 1, 3])

        u_a = np.zeros([N_panels, 3]) # Vetor unitário colinear à corda
        u_n = np.zeros([N_panels, 3]) # Vetor unitário normal à corda

        cp_lengths = np.zeros([N_panels, 3]) # Vetor de comprimento de cada painel
        cp_dsl = np.zeros([N_panels, 3]) # Vetor de comprimento ao longo da envergadura adimensional (dimensionless spanwise length vector)

        cp_areas = np.zeros(N_panels) # Área de cada painel
        cp_chords = np.zeros(N_panels) # Distribuição de cordas dos pontos de colocação
        cp_macs = np.zeros(N_panels) # Corda média aerodinâmica de cada painel
        cp_reynolds = np.zeros(N_panels) # Número de Reynolds de cada painel
        cp_airfoils = np.zeros(N_panels) # Perfil de cada painel

        # TODO: explain what these variables mean
        idx_n = 0
        span_incremental = 0
        height_incremental = 0
        MAC = 0
        for i, span_partition in enumerate(spans):
            n = span_panel_numbers[i]
            chord_i = chords[i]
            chord_ii = chords[i+1]
            offset_i = offsets[i]
            offset_ii = offsets[i+1]
            twist_i = twist_angles[i]
            twist_ii = twist_angles[i+1]
            dihedral_partition = dihedral_angles[i]
            airfoil_i = airfoils[i]
            airfoil_ii = airfoils[i+1]
            # Ângulo de enflechamento [rad] da partição atual
            sweep_partition = np.arctan( (0.25*chord_ii - 0.25*chord_i + offset_i)/span_partition ) if sweep_check is True else 0

            span_y_distr = geo.get_span_y_distr(n, span_partition, distribution_type)
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
                vertice_points[0][0] = vp_x_component[0]
                vertice_points[0][1] = vp_y_component[0]
                vertice_points[0][2] = vp_z_component[0]

            for j, _ in enumerate(cp_y_component):
                # Componente Y dos CP / VP
                vertice_points[idx_n+j+1][1] = span_incremental + vp_y_component[j+1]
                collocation_points[idx_n+j][1] = span_incremental + cp_y_component[j]
                # Componente X
                vertice_points[idx_n+j+1][0] = vp_x_component[j+1]
                collocation_points[idx_n+j][0] = cp_x_component[j]
                # Componente Z
                vertice_points[idx_n+j+1][2] = height_incremental + vp_z_component[j+1]
                collocation_points[idx_n+j][2] = height_incremental + cp_z_component[j]
                # Ângulo de torção (twist geométrico) do ponto de colocação [testar]
                twist_cp = cp_twist_distr[j]

                # Vetores de posicionamento do painel (colinear e perpendicular à corda)
                euler_matrix = geo.get_euler_matrix(dihedral_partition, twist_cp, sweep_partition)
                u_a[idx_n+j] = euler_matrix.dot(np.array([1, 0, 0]))
                u_n[idx_n+j] = euler_matrix.dot(np.array([0, 0, 1]))

                chord_vp_j = geo.get_local_chord(vp_y_component[j], chord_i, chord_ii, span_partition)
                chord_vp_jj = geo.get_local_chord(vp_y_component[j+1], chord_i, chord_ii, span_partition)
                chord_cp = geo.get_local_chord(cp_y_component[j], chord_i, chord_ii, span_partition)
                cp_chords[idx_n+j] = chord_cp

                # Corda média aerodinâmica  do painel
                cp_mac = (2/3)*(chord_vp_j ** 2 + chord_vp_j * chord_vp_jj + chord_vp_jj ** 2)/(chord_vp_j + chord_vp_jj)
                cp_macs[idx_n+j] = cp_mac
                for k in range(3):
                    cp_lengths[idx_n+j][k] = vertice_points[idx_n+j+1][k] - vertice_points[idx_n+j][k]

                # Vetor comprimento do painel
                cp_length_x = cp_lengths[idx_n+j][0]
                cp_length_y = cp_lengths[idx_n+j][1]
                cp_length_z = cp_lengths[idx_n+j][2]

                # Área do painel
                cp_area = 0.5*(chord_vp_j + chord_vp_jj)*math.sqrt( cp_length_y ** 2 + cp_length_z ** 2 )
                cp_areas[idx_n+j] = cp_area
                partition_areas[i] += cp_area

                # Cálculo do dimensionless spanwise length
                for k in range(3): 
                    cp_dsl[idx_n+j][k] = (cp_mac * cp_lengths[idx_n+j][k]) / cp_area

                if airfoil_i == airfoil_ii:
                    # cp_airfoils[idx_n+j] = airfoil_i
                    cp_airfoils[idx_n+j] = 1 # Testar desse jeito depois
                else:
                    cp_airfoils[idx_n+j] = (collocation_points[idx_n+j][1]-vertice_points[idx_n][1]) / span_partition

            # Razão de afilamento da partição atual
            partition_lambda = chord_ii / chord_i
            # Corda média aerodinâmica da asa
            MAC += partition_areas[i]*(2/3)*chord_i*(1+partition_lambda+partition_lambda**2)/(1+partition_lambda)
            idx_n += n
            span_incremental += span_partition
            height_incremental += vp_z_component[-1]
    
        collocation_points[:,0] += self.x_pos
        vertice_points[:,0] += self.x_pos
        collocation_points[:,2] += self.z_pos
        vertice_points[:,2] += self.z_pos
        
        # Atribuição dos valores para o objeto
        self.span_panel_numbers = span_panel_numbers
        self.total_span = total_span
        self.collocation_points = collocation_points
        self.vertice_points = vertice_points
        self.u_a = u_a
        self.u_n = u_n
        self.cp_lengths = cp_lengths
        self.cp_dsl = cp_dsl
        self.cp_areas = cp_areas
        self.cp_chords = cp_chords
        self.cp_macs = cp_macs
        self.cp_reynolds = cp_reynolds
        self.cp_airfoils = cp_airfoils
        self.partition_areas = partition_areas
        self.total_area = sum(partition_areas) * 2
        self.MAC = MAC / sum(partition_areas)
        self.AR = (2 * total_span) ** 2 / (2 * sum(partition_areas))

        def calculate_reynolds(self, flight_condition: FlightCondition):
            pass
