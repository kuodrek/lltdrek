import numpy as np
import math
import models
import matplotlib.pyplot as plt


def get_span_distribution(n, span_partition):
    vertice_points = np.zeros(n + 1)
    collocation_points = np.zeros(n)
    for i in range(n):
        vertice_points[i] = span_partition * 1 / 2 * (1 - np.cos((i * np.pi / n)))
        collocation_points[i] = span_partition * 1 / 2 * (1 - np.cos((i * np.pi / n) - (np.pi / (2 * n))))
    vertice_points[i + 1] = span_partition * 1 / 2 * (1 - np.cos(((i + 1) * np.pi / n)))
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

    # Distribuição dos paineis por partição
    span_panels_distribution = spans / total_span * N_panels
    span_panels_distribution = [math.ceil(int(i)) for i in span_panels_distribution]
    # Adicionar mais um painel por conta dos arredondamentos
    if sum(span_panels_distribution) == N_panels - 1:
        span_panels_distribution[0] += 1

    # Distribuição dos pontos de colocação e vértices
    N_partitions = len(spans)
    
    for i, span_partition in enumerate(spans):
        n = span_panels_distribution[i]
        span_distribution = get_span_distribution(n, span_partition)
        y_cp = np.ones(n)
        y_vp = np.zeros(n+1)
        collocation_points = span_distribution['collocation_points']
        vertice_points = span_distribution['vertice_points']
        plt.scatter(collocation_points, y_cp)
        plt.scatter(vertice_points, y_vp)
        plt.show()
        print(span_distribution)

   

# Teste
b = np.array([3])
asa = models.Wing(
    spans=b,
    chords=0,
    offsets=0,
    twist_angles=0,
    dihedral_angles=0,
    airfoils=0,
    N_panels=20,
    panels_distribution="cosine",
)
generate_mesh(asa)
