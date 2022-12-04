import numpy as np
from models.wing import Wing
from models.flight_condition import FlightCondition
from utils import velocity

asa = Wing(
        spans=[3, 2],
        chords=[1, 1, 1],
        offsets=[0, 0, 0],
        twist_angles=[0, 0, 0],
        dihedral_angles=[0, 0],
        airfoils=['optfoilb2', 'optfoilb2', 'optfoilb2'],
        N_panels=12,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa'
    )

asa.generate_mesh()

asa_2 = Wing(
        spans=[3, 2],
        chords=[1, 1, 1],
        offsets=[0, 0, 0],
        twist_angles=[0, 0, 0],
        dihedral_angles=[0, 0],
        airfoils=['optfoilb2', 'optfoilb2', 'optfoilb2'],
        N_panels=20,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa_2',
        x_pos=1,
        z_pos=1
    )

asa_2.generate_mesh()

flight_condition = FlightCondition(
    V_inf=15,
    nu=1.5e-5,
    rho=1.225,
    aoa = [1],
    ground_effect_check=False
)

v_inf = flight_condition.v_inf_list

# Testing distribution of induced velocities
v_ij_distr_self = velocity.get_induced_velocity_distribution(asa.collocation_points, asa.cp_macs, asa.vertice_points, v_inf)

v_ij_distr = velocity.get_induced_velocity_distribution(asa.collocation_points, asa.cp_macs, asa_2.vertice_points, v_inf)

a = 1