from models.wingpool import WingPool
from models.wing import Wing
from models.flight_condition import FlightCondition


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

flight_condition = FlightCondition(
    V_inf=15,
    nu=1.5e-5,
    rho=1.225,
    aoa = [1],
    ground_effect_check=False
)

wingpool_teste = WingPool([asa], flight_condition)
wingpool_teste.calculate_induced_velocities()
a = 1