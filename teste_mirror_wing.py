import numpy as np
from models.wing import Wing
from models.flight_condition import FlightCondition
from models.wingpool import WingPool
from utils import velocity
import json

# Testing terms of the newton corrector equation and its symmetries

asa = Wing(
        spans=[3],
        chords=[1, 1],
        offsets=[0, 0],
        twist_angles=[0, 0],
        dihedral_angles=[0],
        airfoils=['optfoilb2', 'optfoilb2'],
        N_panels=10,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa',
        x_pos=1,
        z_pos=1
    )

asa_mirrored = Wing(
        spans=[3],
        chords=[1, 1],
        offsets=[0, 0],
        twist_angles=[0, 0],
        dihedral_angles=[0],
        airfoils=['optfoilb2', 'optfoilb2'],
        N_panels=10,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa_mirrored',
        x_pos=0,
        z_pos=0
    )

asa.generate_mesh()
asa_mirrored.generate_mesh()

for idx, cp in enumerate(asa_mirrored.collocation_points):
    asa_mirrored.collocation_points[idx][1] *= -1
    asa_mirrored.vertice_points[idx][1] *= -1
asa_mirrored.vertice_points[idx+1][1] *= -1

flight_condition = FlightCondition(
    V_inf=15,
    nu=1.5e-5,
    rho=1.225,
    aoa = [1],
    ground_effect_check=False
)

wingpool = WingPool(
    wing_list=[asa, asa_mirrored],
    flight_condition=flight_condition
)

induced_velocities = wingpool.calculate_induced_velocities()
total_velocity = wingpool.get_total_dim_velocity(asa)


## Inverter 
with open("teste_simetria_asas_2.txt", "w") as f:
    asa_asa = induced_velocities['asa']['asa'][0]
    asa_asa_mirrored = induced_velocities['asa']['asa_mirrored'][0]
    asa_mirrored_asa_mirrored = induced_velocities['asa_mirrored']['asa_mirrored'][0]
    asa_mirrored_asa = induced_velocities['asa_mirrored']['asa'][0]
    f.write("\nasa - asa\n")
    for velocity in asa_asa:
        f.write(f"{velocity[0]},{velocity[1]},{velocity[2]}\n")

    f.write("\nasa - asa mirrored\n")
    for velocity in asa_asa_mirrored:
        f.write(f"{velocity[0]},{velocity[1]},{velocity[2]}\n")

    f.write("\nasa mirrored - asa mirrored\n")
    for velocity in asa_mirrored_asa_mirrored:
        f.write(f"{velocity[0]},{velocity[1]},{velocity[2]}\n")
    
    f.write("\nasa mirrored - asa\n")
    for velocity in asa_mirrored_asa:
        f.write(f"{velocity[0]},{velocity[1]},{velocity[2]}\n")