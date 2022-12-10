import numpy as np
from models.wing import Wing
from models.flight_condition import FlightCondition
from models.wingpool import WingPool
from utils import velocity

# Testing terms of the newton corrector equation and its symmetries

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
flight_condition = FlightCondition(
    V_inf=15,
    nu=1.5e-5,
    rho=1.225,
    aoa = [1],
    ground_effect_check=False
)

wingpool = WingPool(
    wing_list=[asa],
    flight_condition=flight_condition
)

induced_velocities = wingpool.calculate_induced_velocities()
total_velocity = wingpool.get_total_dim_velocity(asa)

# Original wing
u_n = asa.u_n.copy()
u_a = asa.u_a.copy()
cp_dsl = asa.cp_dsl.copy()

w_list = []
for idx, dsl in enumerate(cp_dsl):
    w_list.append(np.cross(total_velocity, dsl))

# Symmetry wing
w_list_sym = []
for idx, dsl in enumerate(cp_dsl):
    dsl[1] = dsl[1]*-1

    w_list_sym.append(np.cross(total_velocity, dsl))

for i, _ in enumerate(w_list):
    print("------------------------")
    print(f"original: {w_list[i]}")
    print(f"symmetric: {w_list_sym[i]}")