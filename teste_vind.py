import numpy as np
from models.wing import Wing
# from utils.velocity import get_induced_velocity

asa = Wing(
        spans=[3, 2],
        chords=[1, 0.8, 0.4],
        offsets=[0, 0, 0.5],
        twist_angles=[0, 0, 0],
        dihedral_angles=[10, 15],
        airfoils=['optfoilb2', 'optfoilb2', 'optfoilb2'],
        N_panels=12,
        distribution_type="cosine",
        sweep_check=False
    )

print(f'antes de generate mesh: {asa.u_a}')
asa.generate_mesh()
print(f'depois de generate mesh: {asa.u_a}')
# collocation_point, vertice_point_1, vertice_point_2, v_inf, mac, h, same_panel_check)
# cp = asa.collocation_points[0]
# vp_1 = asa.vertice_points[0]
# vp_2 = asa.vertice_points[1]
# MAC = asa.MAC
# v_inf = np.array([0.9947, 0, -0.1028])
# get_induced_velocity(cp, vp_1, vp_2, v_inf, MAC)