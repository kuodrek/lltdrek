from xml.etree.ElementTree import C14NWriterTarget
import numpy as np
import math

from pyparsing import cpp_style_comment
import models
import matplotlib.pyplot as plt


def plot_points(collocation_points, vertice_points):
    figure, axis = plt.subplots(2, 1)
    axis[0].invert_yaxis()
    axis[0].scatter(vertice_points[:,1], vertice_points[:,0], label='vertices')
    axis[0].scatter(collocation_points[:,1], collocation_points[:,0], label='collocation points')
    axis[0].set_title("X x Y")
    axis[0].legend(loc='best')

    # For Cosine Function
    axis[1].scatter(vertice_points[:,1], vertice_points[:,2], label='vertices')
    axis[1].scatter(collocation_points[:,1], collocation_points[:,2], label='collocation points')
    axis[1].set_title("Z x Y")
    axis[1].legend(loc='best')
    plt.show()


def get_span_y_distr(n, span_partition, distribution_type):
    vertice_points = np.zeros(n + 1)
    collocation_points = np.zeros(n)

    vertice_points[0] = span_partition * 1 / 2 * \
        (1 - np.cos(0))  # Note that this is = 0
    for i in range(n):
        # Y component of points
        if distribution_type == 'cosine':
            vertice_points[i+1] = span_partition * 1 / \
                2 * (1 - np.cos(((i+1) * np.pi / n)))
            collocation_points[i] = span_partition * 1 / 2 * \
                (1 - np.cos(((i+1) * np.pi / n) - (np.pi / (2 * n))))
        elif distribution_type == 'linear':
            vertice_points[i+1] = span_partition * ((i+1)/n)
            collocation_points[i] = span_partition * ((i)/n + (i+1)/n) * 1 / 2
        else:
            raise Exception(
                "Invalid 'distribution_type' value. Expected 'cosine' or 'linear'.")

    span_y_distribution = {
        'vertice_points': vertice_points,
        'collocation_points': collocation_points
    }
    return span_y_distribution


def get_span_x_distr(cp_y_component, vp_y_component, chord_i, chord_ii, offset_i, offset_ii, span_partition):
    n = len(cp_y_component)
    vp_x_component = np.zeros(n+1)
    cp_x_component = np.zeros(n)

    vp_x_component[0] = \
    0.25*(chord_ii - chord_i) / span_partition * vp_y_component[0] + 0.25 * chord_i + \
    (offset_ii - offset_i) / span_partition * vp_y_component[0] + offset_i
    for i in range(n):
        vp_x_component[i+1] = \
        0.25*(chord_ii - chord_i) / span_partition * vp_y_component[i+1] + 0.25 * chord_i + \
        (offset_ii - offset_i) / span_partition * vp_y_component[i+1] + offset_i

        cp_x_component[i] = \
        0.25*(chord_ii - chord_i) / span_partition * cp_y_component[i] + 0.25 * chord_i + \
        (offset_ii - offset_i) / span_partition * cp_y_component[i] + offset_i
    
    span_x_distribution = {
        'vertice_points': vp_x_component,
        'collocation_points': cp_x_component
    }
    return span_x_distribution


def get_span_z_distr(cp_y_component, vp_y_component, dihedral_partition):
    n = len(cp_y_component)
    vp_z_component = np.zeros(n+1)
    cp_z_component = np.zeros(n) 

    vp_z_component[0] = np.tan(dihedral_partition) * vp_y_component[0]
    for i in range(n):
        vp_z_component[i+1] = np.tan(dihedral_partition) * vp_y_component[i+1]
        cp_z_component[i] = np.tan(dihedral_partition) * cp_y_component[i]
    
    span_z_distribution = {
        'vertice_points': vp_z_component,
        'collocation_points': cp_z_component
    }
    return span_z_distribution


def get_twist_distr(cp_y_component, twist_i, twist_ii, span_partition):
    n = len(cp_y_component)
    span_twist_distribution = np.zeros(n)
    for i in range(n):
        span_twist_distribution[i] = (twist_ii - twist_i) / span_partition * cp_y_component[i] + twist_i
    return span_twist_distribution


def get_euler_matrix(dihedral, twist, sweep):
    # Utilizando-se rotação ZYX -> Ângulos de Tait-Bryan
    # REFERÊNCIA: https://en.wikipedia.org/wiki/Euler_angles#Conversion_to_other_orientation_representations
    c1 = np.cos(dihedral)
    c2 = np.cos(twist)
    c3 = np.cos(sweep)
    s1 = np.sin(dihedral)
    s2 = np.sin(twist)
    s3 = np.sin(sweep)
    euler_matrix = np.array([
        [c1*c2, c1*s2*s3-c3*s1, s1*s3+c1*c3*s2],
        [c2*s1, c1*c3+s1*s2*s3, c3*s1*s2-c1*s3],
        [-s2, c2*s3, c2*c3]
        ])
    return euler_matrix


def generate_mesh(Wing: models.Wing):
    # Inicialização das variáveis
    N_panels = Wing.N_panels
    spans = Wing.spans
    chords = Wing.chords
    offsets = Wing.offsets
    twist_angles = Wing.twist_angles
    dihedral_angles = Wing.dihedral_angles
    distribution_type = Wing.distribution_type

    total_span = sum(spans)
    # Número de painéis por partição
    span_panel_numbers = spans / total_span * N_panels
    span_panel_numbers = [math.ceil(int(i)) for i in span_panel_numbers]
    # Adicionar mais um painel por conta dos arredondamentos
    if sum(span_panel_numbers) == N_panels - 1:
        span_panel_numbers[-1] += 1

    # Distribuição dos pontos de colocação e vértices
    # Note que o ponto de colocação é a representação pontual de um painel
    collocation_points = np.zeros([N_panels, 3])
    vertice_points = np.zeros([N_panels + 1, 3])

    # Distribuição de vetores unitários dos pontos de colocação
    u_a = np.zeros([N_panels, 3]) # Vetor unitário no sentido da corda
    u_n = np.zeros([N_panels, 3]) # Vetor unitário normal à corda

    # Distribuição de propriedades geométricas dos painéis da asa
    cp_lengths = np.zeros([N_panels, 3]) # Vetor de comprimento de cada painel
    cp_dsl = np.zeros([N_panels, 3]) # Vetor de comprimento ao longo da envergadura adimensional (dimensionless spanwise length vector)

    # Distribuição de informações geométricas gerais dos pontos de colocação
    cp_areas = np.zeros(N_panels) # Área de cada painel
    cp_mac = np.zeros(N_panels) # Corda média aerodinâmica de cada painel
    cp_reynolds = np.zeros(N_panels) # Número de Reynolds de cada painel
    cp_airfoil = np.zeros(N_panels) # Perfil de cada painel

    # Explicar estas variáveis
    idx_i = 0
    span_incremental = 0
    height_incremental = 0
    for i, span_partition in enumerate(spans):
        n = span_panel_numbers[i]
        chord_i = chords[i]
        chord_ii = chords[i+1]
        offset_i = offsets[i]
        offset_ii = offsets[i+1]
        twist_i = twist_angles[i]
        twist_ii = twist_angles[i+1]
        dihedral_partition = dihedral_angles[i]
        # Ângulo de enflechamento [rad] da partição atual
        sweep_partition = np.arctan( (0.25*chord_ii - 0.25*chord_i + offset_i)/span_partition )

        span_y_distr = get_span_y_distr(n, span_partition, distribution_type)
        cp_y_component = span_y_distr['collocation_points']
        vp_y_component = span_y_distr['vertice_points']

        span_x_distr = get_span_x_distr(cp_y_component,vp_y_component,chord_i,chord_ii,offset_i,offset_ii,span_partition)
        cp_x_component = span_x_distr['collocation_points']
        vp_x_component = span_x_distr['vertice_points']

        span_z_distr = get_span_z_distr(cp_y_component,vp_y_component,dihedral_partition)
        cp_z_component = span_z_distr['collocation_points']
        vp_z_component = span_z_distr['vertice_points']

        cp_twist_distr = get_twist_distr(cp_y_component, twist_i, twist_ii, span_partition)

        for j, _ in enumerate(cp_y_component):
            # Componente Y dos CP / VP
            collocation_points[idx_i+j][1] = span_incremental + cp_y_component[j]
            vertice_points[idx_i+j][1] = span_incremental + vp_y_component[j]
            # Componente X
            collocation_points[idx_i+j][0] = cp_x_component[j]
            vertice_points[idx_i+j][0] = vp_x_component[j]
            # Componente Z
            collocation_points[idx_i+j][2] = height_incremental + cp_z_component[j]
            vertice_points[idx_i+j][2] = height_incremental + vp_z_component[j]

            # Ângulo de torção (twist geométrico) do ponto de colocação
            twist_cp = cp_twist_distr[j]
            euler_matrix = get_euler_matrix(dihedral_partition, twist_cp, sweep_partition)
            u_a[idx_i+j] = euler_matrix.dot(np.array([1, 0, 0]))
            u_n[idx_i+j] = euler_matrix.dot(np.array([0, 0, 1]))

            # Implementar:
            # cp_lengths
            # cp_dsl
            
            # cp_areas
            # cp_mac
            # cp_reynolds
            # cp_airfoil
        vertice_points[idx_i+j+1][1] = span_incremental + vp_y_component[-1]
        vertice_points[idx_i+j+1][0] = vp_x_component[-1]
        vertice_points[idx_i+j+1][2] = height_incremental + vp_z_component[-1]

        idx_i += n
        span_incremental += span_partition
        height_incremental += vp_z_component[-1]




# Teste
asa = models.Wing(
    spans=[3, 2],
    chords=[1, 0.8, 0.4],
    offsets=[0, 0, 0.5],
    twist_angles=[0, 0, 0],
    dihedral_angles=[10, 15],
    airfoils=['optfoilb2', 'optfoilb2', 'optfoilb2'],
    N_panels=10,
    distribution_type="cosine",
)
generate_mesh(asa)
