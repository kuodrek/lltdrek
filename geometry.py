import numpy as np
import math
import models
import matplotlib.pyplot as plt


def get_span_distribution(n, span_partition, distribution_type):
    vertice_points = np.zeros(n + 1)
    collocation_points = np.zeros(n)

    vertice_points[0] = span_partition * 1 / 2 * (1 - np.cos(0)) # Note that this is = 0
    for i in range(n):
        # Y component of points
        if distribution_type == 'cosine':
            vertice_points[i+1] = span_partition * 1 / 2 * (1 - np.cos(((i+1) * np.pi / n)))
            collocation_points[i] = span_partition * 1 / 2 * (1 - np.cos(((i+1) * np.pi / n) - (np.pi / (2 * n))))
        elif distribution_type == 'linear':
            vertice_points[i+1] = span_partition * ((i+1)/n)
            collocation_points[i] = span_partition * ((i)/n + (i+1)/n) * 1 / 2         
        else: raise Exception("Invalid 'distribution_type' value. Expected 'cosine' or 'linear'.")
    
    span_distribution = {
        'vertice_points': vertice_points,
        'collocation_points': collocation_points
    }
    return span_distribution


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
    collocation_points = np.zeros([N_panels,3])
    vertice_points = np.zeros([N_panels+1,3])

    # Componente no eixo X dos pontos


    # Componente no eixo Y dos pontos
    idx_i = 0
    span_incremental = 0
    # Vetores para verificação, remover depois
    vp_teste = []
    cp_teste = []
    for i, span_partition in enumerate(spans):
        n = span_panels_distribution[i]

        span_distribution = get_span_distribution(n, span_partition, distribution_type)
        collocation_points_partition = span_distribution['collocation_points']
        vertice_points_partition = span_distribution['vertice_points']

        cp_teste.append(collocation_points_partition)
        vp_teste.append(vertice_points_partition)

        for j, _ in enumerate(collocation_points_partition):
            collocation_points[idx_i+j][1] = span_incremental + collocation_points_partition[j]
            vertice_points[idx_i+j][1] = span_incremental + vertice_points_partition[j]
        vertice_points[idx_i+j+1][1] = span_incremental + vertice_points_partition[-1]

        idx_i += n
        span_incremental += spans[i]

   

# Teste
b = np.array([3,2,1])
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
