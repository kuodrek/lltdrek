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

wingpool.calculate_induced_velocities()
wingpool.calculate_aoa_eff()

a=1