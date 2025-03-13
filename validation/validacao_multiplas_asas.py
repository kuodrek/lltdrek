from testes.validacao_aoas import aoa_list_1
from lltdrek.models.wing import Wing
from lltdrek.models.flight_condition import FlightCondition
from lltdrek.models.wingpool import WingPool
from lltdrek.models.simulation import Simulation
from lltdrek.models.post_processing import PostProcessing
from lltdrek.utils import data


airfoils_data_dict, airfoils_dat_dict = data.load_folder('airfoils_test')

# Testing terms of the newton corrector equation and its symmetries
# print("ATENÇAO, AS FUNCOES DE LOOKUP ESTAO RETORNANDO 1 PARA CL E CL_ALPHA. CL0, CM0 PARA FINS DE COMPARACAO")
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

asa_2 = Wing(
        spans=[4.572/2],
        chords=[0.813, 0.325],
        offsets=[0, 0.122],
        twist_angles=[0, -2.4],
        dihedral_angles=[0, 0],
        airfoils=['NACA4412_exp', 'NACA4424_exp'],
        N_panels=12,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa_2',
        x_pos=0,
        z_pos=0.5
    )

asa.generate_mesh()
asa_2.generate_mesh()

flight_condition = FlightCondition(
    V_inf=45.48,
    nu=6.39e-6,
    rho=2.83,
    # aoa = [aoa_list_1[0]],
    # aoa = aoa_list_1,
    aoa = [13, 14, 15, 16, 17, 18, 19, 20, 21],
    ground_effect_check=False
)

asa.setup_airfoil_data(flight_condition, airfoils_data_dict)
asa_2.setup_airfoil_data(flight_condition, airfoils_data_dict)

wingpool = WingPool(
    wing_list=[asa, asa_2],
    flight_condition=flight_condition
)

simulation = Simulation(
    wing_pool=wingpool,
    damping_factor=1,
    max_iter=150,
    max_residual=1e-3,
    linear_check=True
)

G_solution_list = simulation.run_simulation()

llt_coeficientes = PostProcessing(ref_point=[0, 0, 0])

coefs_list = []
coefs_single_list = []
for i, G in enumerate(G_solution_list):
    coefs = llt_coeficientes.get_global_coefficients(wingpool, G, aoa_index=i, S_ref=wingpool.S_ref, c_ref=asa.MAC)
    coefs_single = llt_coeficientes.get_wing_coefficients(wingpool, G, aoa_index=i, S_ref=asa.total_area*2, c_ref=asa.MAC)

    coefs_list.append(coefs["CF"][2])
    coefs_single_list.append([coefs_single["asa"]["CF"][2], coefs_single["asa_2"]["CF"][2]])

CLmax_dict = llt_coeficientes.get_CL_max_linear(wing_pool=wingpool, G_list=G_solution_list, S_ref=wingpool.S_ref)
a=1