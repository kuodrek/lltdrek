import numpy as np
from models.wing import Wing
from models.flight_condition import FlightCondition
from utils import velocity
import copy

asa = Wing(
        spans=[3, 2],
        chords=[1, 0.8, 0.4],
        offsets=[0, 0, 0.5],
        twist_angles=[0, 0, 0],
        dihedral_angles=[10, 15],
        airfoils=['optfoilb2', 'optfoilb2', 'optfoilb2'],
        N_panels=12,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa'
    )

asa.generate_mesh()
asa = copy.deepcopy(asa)
print(len(asa.collocation_points))

asa_2 = Wing(
        spans=[3, 2],
        chords=[1, 0.8, 0.4],
        offsets=[0, 0, 0.5],
        twist_angles=[0, 0, 0],
        dihedral_angles=[10, 15],
        airfoils=['optfoilb2', 'optfoilb2', 'optfoilb2'],
        N_panels=20,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa_2',
        x_pos=1,
        z_pos=1
    )

asa_2.generate_mesh()
print(len(asa_2.collocation_points))

flight_condition = FlightCondition(
    V_inf=15,
    nu=1.5e-5,
    rho=1.225,
    aoa = [1],
    ground_effect_check=False
)

cp = asa.collocation_points[0]
vp_1 = asa.vertice_points[0]
vp_2 = asa.vertice_points[1]
cp_mac = asa.cp_macs[0]
v_inf = flight_condition.v_inf_array

# Testing single induced velocity
velocity.get_induced_velocity(cp, vp_1, vp_2, cp_mac, v_inf)

# Testing distribution of induced velocities
v_ij_distr = velocity.get_induced_velocity_distribution(asa.collocation_points, asa.cp_macs, asa_2.vertice_points, v_inf)

a = 1