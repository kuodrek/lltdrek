from validacao_aoas import aoa_list_2
from models.wing import Wing
from models.flight_condition import FlightCondition
from models.wingpool import WingPool
from models.simulation import Simulation
from utils import data


airfoils_data_dict, airfoils_dat_dict = data.load_folder('airfoils_test')

asa = Wing(
        spans=[4.572/2],
        chords=[0.726, 0.290],
        offsets=[0, 0.109],
        twist_angles=[0, -2],
        dihedral_angles=[0, 3],
        airfoils=['NACA64210_xfoil', 'NACA64210_xfoil'],
        N_panels=12,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa'
    )

asa.generate_mesh()

flight_condition = FlightCondition(
    V_inf=57.85,
    nu=7.44e-6,
    rho=2.43,
    aoa = aoa_list_2,
    ground_effect_check=False
)

asa.setup_airfoil_data(flight_condition, airfoils_data_dict)

wingpool = WingPool(
    wing_list=[asa],
    flight_condition=flight_condition
)

simulation = Simulation(
    wing_pool=wingpool,
    damping_factor=0.8,
    max_iter=1000,
    max_residual=1e-3,
    linear_check=False
)

G_solution_list = simulation.run_simulation()

a=1