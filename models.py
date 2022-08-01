import numpy as np
import geometry as geo
import math
"""
inputs iniciais
ASA
Vetor de envergaduras - b
vetor de cordas - c
vetor de offsets - offsets
vetor de twists - twangles
vetor de diedros - dihedral
vetor de perfis - airfoils
gap - para o caso de duas asas

CONDICOES DE VOO
Velocidade do escoamento - Vinf
Coeficiente de viscosidade dinamica do ar - nu
densidade do ar - rho
altura do chao - h
angulo de ataque - aoa
Numero de paineis - N
Distribuicao dos paineis - panel_distribution

SIMULACAO
Fator de amortecimento - damp
Fator de convergencia - R
Numero de iteracoes - max_iter
"""


class Wing:
    def __init__(
        self,
        spans,
        chords,
        offsets,
        twist_angles,
        dihedral_angles,
        airfoils,
        N_panels=20,
        distribution_type="linear",
        sweep_check=False
    ):
        Wing.spans = spans
        Wing.chords = chords
        Wing.offsets = offsets
        Wing.twist_angles = [angle * np.pi / 180 for angle in twist_angles]
        Wing.dihedral_angles = [angle * np.pi / 180 for angle in dihedral_angles]
        Wing.airfoils = airfoils
        Wing.N_panels = N_panels
        Wing.distribution_type = distribution_type
        Wing.sweep_check = sweep_check
        # Propriedades obtidas através do método get_mesh
        Wing.total_span = None
        Wing.span_panel_numbers = None
        Wing.collocation_points = None
        Wing.vertice_points = None
        Wing.u_a = None
        Wing.u_n = None
        Wing.cp_lengths = None
        Wing.cp_dsl = None
        Wing.cp_areas = None
        Wing.cp_chords = None
        Wing.cp_macs = None
        Wing.cp_reynolds = None
        Wing.cp_airfoils = None
        Wing.MAC = None
        Wing.partition_areas = None
        Wing.total_area = None
        Wing.AR = None
    
    def generate_mesh(self):
        # Inicialização das variáveis
        N_panels = self.N_panels
        spans = self.spans
        chords = self.chords
        offsets = self.offsets
        twist_angles = self.twist_angles
        dihedral_angles = self.dihedral_angles
        airfoils = self.airfoils
        distribution_type = self.distribution_type
        sweep_check = self.sweep_check

        self.total_span = sum(spans)
        # Número de painéis por partição
        # span_panel_numbers = spans / self.total_span * N_panels
        self.span_panel_numbers = [value / self.total_span * N_panels for value in spans]
        self.span_panel_numbers = [math.ceil(int(i)) for i in self.span_panel_numbers]
        # Adicionar mais um painel por conta dos arredondamentos
        if sum(self.span_panel_numbers) == N_panels - 1:
            self.span_panel_numbers[-1] += 1

        self.partition_areas = np.zeros(len(spans))
        # Distribuição dos pontos de colocação e vértices
        # Note que o ponto de colocação é a representação pontual de um painel
        self.collocation_points = np.zeros([N_panels, 3])
        self.vertice_points = np.zeros([N_panels + 1, 3])

        # Distribuição de vetores unitários dos pontos de colocação (posição da seção em relação ao escoamento)
        self.u_a = np.zeros([N_panels, 3]) # Vetor unitário colinear à corda
        self.u_n = np.zeros([N_panels, 3]) # Vetor unitário normal à corda

        # Distribuição de propriedades geométricas dos painéis da asa
        self.cp_lengths = np.zeros([N_panels, 3]) # Vetor de comprimento de cada painel
        self.cp_dsl = np.zeros([N_panels, 3]) # Vetor de comprimento ao longo da envergadura adimensional (dimensionless spanwise length vector)

        # Distribuição de informações geométricas gerais dos pontos de colocação
        self.cp_areas = np.zeros(N_panels) # Área de cada painel
        self.cp_chords = np.zeros(N_panels) # Distribuição de cordas dos pontos de colocação
        self.cp_macs = np.zeros(N_panels) # Corda média aerodinâmica de cada painel
        self.cp_reynolds = np.zeros(N_panels) # Número de Reynolds de cada painel
        self.cp_airfoils = np.zeros(N_panels) # Perfil de cada painel

        # Explicar estas variáveis
        idx_n = 0
        span_incremental = 0
        height_incremental = 0
        MAC = 0
        for i, span_partition in enumerate(spans):
            n = self.span_panel_numbers[i]
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
                self.vertice_points[0][1] = vp_y_component[0]
                self.vertice_points[0][0] = vp_x_component[0]
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

                # cálculo do dimensionless spanwise length
                for k in range(3): 
                    self.cp_dsl[idx_n+j][k] = (cp_mac * self.cp_lengths[idx_n+j][k]) / cp_area

                if airfoil_i == airfoil_ii:
                    # self.cp_airfoils[idx_n+j] = airfoil_i
                    self.cp_airfoils[idx_n+j] = 1 # Testar desse jeito depois
                else:
                    self.cp_airfoils[idx_n+j] = (self.collocation_points[idx_n+j][1]-self.vertice_points[idx_n][1]) / span_partition

            # Razão de afilamento da partição atual
            partition_lambda = chord_ii / chord_i
            # Corda média aerodinâmica da asa
            MAC += self.partition_areas[i]*(2/3)*chord_i*(1+partition_lambda+partition_lambda**2)/(1+partition_lambda)
            idx_n += n
            span_incremental += span_partition
            height_incremental += vp_z_component[-1]
    
        self.MAC = MAC / sum(self.partition_areas)
        self.total_area = sum(self.partition_areas) * 2
        self.AR = (2 * self.total_span) ** 2 / (2 * sum(self.partition_areas))

class FlightConditions:
    def __init__(self, V_inf, nu, rho, h=0):
        FlightConditions.V_inf = V_inf
        FlightConditions.nu = nu
        FlightConditions.rho = rho
        FlightConditions.h = h


class Simulation:
    def __init__(self, damp=0.7, max_iter=1000, residual_value=0.001):
        Simulation.damp = damp
        Simulation.max_iter = max_iter
        Simulation.residual_value = residual_value
