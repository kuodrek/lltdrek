import numpy as np
import math
import models
import matplotlib.pyplot as plt


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


def get_span_x_distr(n, cp_y_component, vp_y_component, chord_i, chord_f, offset_i, offset_f, span_partition):
    cp_x_component = np.zeros(n+1)
    vp_x_component = np.zeros(n)

    vp_x_component[0] = \
    0.25*(chord_f - chord_i)/(span_partition) * (vp_y_component[0]-0.25*chord_i) + 0.25*chord_i + \
    (offset_f - offset_i)/(span_partition) * (vp_y_component[0]-offset_i) + offset_i
    for i in range(n):
        vp_x_component[i+1] = \
        0.25*(chord_f - chord_i)/(span_partition) * (vp_y_component[i+1]-0.25*chord_i) + 0.25*chord_i + \
        (offset_f - offset_i)/(span_partition) * (vp_y_component[i+1]-offset_i) + offset_i

        cp_x_component[i] = \
        0.25*(chord_f - chord_i)/(span_partition) * (cp_y_component[i]-0.25*chord_i) + 0.25*chord_i + \
        (offset_f - offset_i)/(span_partition) * (cp_y_component[i]-offset_i) + offset_i
    
    span_x_distribution = {
        'vertice_points': vp_x_component,
        'collocation_points': cp_x_component
    }
    return span_x_distribution


def get_span_z_distr():
    pass


def generate_mesh(Wing: models.Wing):
    # Inicialização das variáveis
    N_panels = Wing.N_panels
    spans = Wing.spans
    total_span = sum(spans)
    chords = Wing.chords
    offsets = Wing.offsets
    dihedral_angles = Wing.dihedral_angles
    distribution_type = Wing.distribution_type

    # Distribuição dos paineis por partição
    span_panels_distribution = spans / total_span * N_panels
    span_panels_distribution = [math.ceil(int(i)) for i in span_panels_distribution]
    # Adicionar mais um painel por conta dos arredondamentos
    if sum(span_panels_distribution) == N_panels - 1:
        span_panels_distribution[0] += 1

    # Distribuição dos pontos de colocação e vértices
    collocation_points = np.zeros([N_panels, 3])
    vertice_points = np.zeros([N_panels+1, 3])

    # Componente no eixo X dos pontos

    idx_i = 0
    span_incremental = 0
    offset_incremental = 0
    height_incremental = 0
    # Vetores para verificação, remover depois
    for i, span_partition in enumerate(spans):
        n = span_panels_distribution[i]
        chord_i = 1
        chord_f = 2
        offset_i = 1
        offset_f = 2

        span_y_distr = get_span_y_distr(n, span_partition, distribution_type)
        cp_y_component = span_y_distr['collocation_points']
        vp_y_component = span_y_distr['vertice_points']

        span_x_distr = get_span_x_distr(n,cp_y_component,vp_y_component,chord_i,chord_f,offset_i,offset_f,span_partition)
        cp_x_component = span_x_distr['collocation_points']
        vp_x_component = span_x_distr['vertice_points']

        # span_z_distr = get_span_z_distr(n,cp_y_component,vp_y_component,dihedral_i,dihedral_f,span_partition)
        # cp_z_component = span_z_distr['collocation_points']
        # vp_z_component = span_z_distr['vertice_points']

        for j, _ in enumerate(cp_y_component):
            # Componente Y dos CP / VP
            collocation_points[idx_i+j][1] = span_incremental + cp_y_component[j]
            vertice_points[idx_i+j][1] = span_incremental + vp_y_component[j]
            # Componente X
            collocation_points[idx_i+j][0] = offset_incremental + cp_x_component[j]
            vertice_points[idx_i+j][0] = offset_incremental + vp_x_component[j]
            # Componente Z
        vertice_points[idx_i+j+1][1] = span_incremental + vp_y_component[-1]
        vertice_points[idx_i+j+1][0] = offset_incremental + vp_x_component[-1]

        idx_i += n
        span_incremental += spans[i]
        offset_incremental += offsets[i]


# Teste
b = np.array([3, 2, 1])
asa = models.Wing(
    spans=b,
    chords=0,
    offsets=0,
    twist_angles=0,
    dihedral_angles=0,
    airfoils=0,
    N_panels=13,
    distribution_type="cosine",
)
generate_mesh(asa)
