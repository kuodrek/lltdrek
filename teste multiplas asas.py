from validacao_aoas import aoa_list_1
from models.wing import Wing
from models.flight_condition import FlightCondition
from models.wingpool import WingPool
from models.simulation import Simulation
from utils import data


airfoils_data_dict, airfoils_dat_dict = data.load_folder('airfoils_test')

# Testing terms of the newton corrector equation and its symmetries

asa_baixo = Wing(
        spans=[4.572/2],
        chords=[0.813, 0.325],
        offsets=[0, 0.122],
        twist_angles=[0, -2.4],
        dihedral_angles=[0, 0],
        airfoils=['NACA4412_exp', 'NACA4424_exp'],
        N_panels=12,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa_baixo'
    )

asa_cima = Wing(
        spans=[4.572/2],
        chords=[0.813, 0.325],
        offsets=[0, 0.122],
        twist_angles=[0, -2.4],
        dihedral_angles=[0, 0],
        airfoils=['NACA4412_exp', 'NACA4424_exp'],
        N_panels=12,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa_cima',
        x_pos=3,
        z_pos=3,
    )

asa_baixo.generate_mesh()
asa_cima.generate_mesh()

flight_condition = FlightCondition(
    V_inf=45.48,
    nu=6.39e-6,
    rho=2.83,
    aoa = aoa_list_1,
    # aoa=[1, 2, 3],
    ground_effect_check=False
)

asa_baixo.setup_airfoil_data(flight_condition, airfoils_data_dict)
asa_cima.setup_airfoil_data(flight_condition, airfoils_data_dict)

wingpool = WingPool(
    wing_list=[asa_baixo, asa_cima],
    flight_condition=flight_condition
)

simulation = Simulation(
    wing_pool=wingpool,
    damping_factor=0.05,
    max_iter=1000,
    max_residual=1e-3,
    linear_check=False
)

G_solution_list = simulation.run_simulation()

a=1