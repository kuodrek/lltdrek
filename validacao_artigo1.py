from validacao_aoas import aoa_list_1
from models.wing import Wing
from models.flight_condition import FlightCondition
from models.wingpool import WingPool
from models.simulation import Simulation
from models.post_processing import PostProcessing
from utils import data


airfoils_data_dict, airfoils_dat_dict = data.load_folder('airfoils_test')

# Testing terms of the newton corrector equation and its symmetries
print("ATENÇAO, AS FUNCOES DE LOOKUP ESTAO RETORNANDO 1 PARA CL E CL_ALPHA. CL0, CM0 PARA FINS DE COMPARACAO")
asa = Wing(
        spans=[4.572/2],
        chords=[0.813, 0.325],
        offsets=[0, 0.122],
        twist_angles=[0, -2.4],
        dihedral_angles=[0, 0],
        airfoils=['NACA4412_exp', 'NACA4424_exp'],
        N_panels=12,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa'
    )


asa.generate_mesh()

flight_condition = FlightCondition(
    V_inf=45.48,
    nu=6.39e-6,
    rho=2.83,
    aoa = [aoa_list_1[0]],
    # aoa = aoa_list_1,
    ground_effect_check=False
)

asa.setup_airfoil_data(flight_condition, airfoils_data_dict)

wingpool = WingPool(
    wing_list=[asa],
    flight_condition=flight_condition
)

simulation = Simulation(
    wing_pool=wingpool,
    damping_factor=0.2,
    max_iter=1000,
    max_residual=1e-3,
    linear_check=False
)

G_solution_list = simulation.run_simulation()

llt_coeficientes = PostProcessing(ref_point=[0, 0, 0])

# coefs = llt_coeficientes.get_global_coefficients(wingpool, G_solution_list[0], aoa_index=0, S_ref=asa.total_area*2, c_ref=asa.MAC)
CLmax = llt_coeficientes.get_CL_max_linear(wingpool, G_solution_list, None, S_ref=asa.total_area*2)
a=1