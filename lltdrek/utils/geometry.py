import numpy as np
import matplotlib.pyplot as plt
from lltdrek.models.enums import DistributionTypes


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


def get_span_y_distr(n: int, span_partition: int, distribution_type: DistributionTypes) -> dict:
    vertice_points = np.zeros(n + 1)
    collocation_points = np.zeros(n)

    vertice_points[0] = span_partition * 1 / 2 * \
        (1 - np.cos(0))  # O primeiro ponto é 0
    for i in range(n):
        # Componente Y dos pontos
        if distribution_type == DistributionTypes.Cosine:
            vertice_points[i+1] = span_partition * 1 / \
                2 * (1 - np.cos(((i+1) * np.pi / n)))
            collocation_points[i] = span_partition * 1 / 2 * \
                (1 - np.cos(((i+1) * np.pi / n) - (np.pi / (2 * n))))
        elif distribution_type == DistributionTypes.Linear:
            vertice_points[i+1] = span_partition * ((i+1)/n)
            collocation_points[i] = span_partition * ((i)/n + (i+1)/n) * 1 / 2
        else:
            raise Exception(
                f"Invalid distribution type. Allowed values: {print(DistributionTypes)}"
            )

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


def get_local_chord(y, chord_i, chord_ii, span_partition):
    return (chord_ii - chord_i) / span_partition * y + chord_i

def get_euler_matrix(dihedral, twist, sweep):
    # Utilizando-se rotação XYZ -> Ângulos de Tait-Bryan
    # REFERÊNCIA: https://en.wikipedia.org/wiki/Euler_angles#Conversion_to_other_orientation_representations
    c1 = np.cos(dihedral)
    c2 = np.cos(twist)
    c3 = np.cos(sweep)
    s1 = np.sin(dihedral)
    s2 = np.sin(twist)
    s3 = np.sin(sweep)
    euler_matrix = np.array([
        [c2*c3, -c2*s3, s2],
        [c1*s3+c3*s1*s2, c1*c3-s1*s2*s3, -c2*s1],
        [s1*s3-c1*c3*s2, c3*s1+c1*s2*s3, c1*c2]
        ])
    return euler_matrix
