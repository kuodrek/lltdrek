from lltdrek.models.wing import Wing
from lltdrek.models.flight_condition import FlightCondition
from lltdrek.models.wingpool import WingPool
from lltdrek.models.simulation import Simulation
from lltdrek.models.post_processing import PostProcessing
from lltdrek.utils import data
import numpy as np

airfoils_data_dict, airfoils_dat_dict = data.load_folder('perfis_clmax_bugado', 0, 8)

# Testando eh bugando na micro
# marina03_alterado: cortado reynolds que ta com dados cagados
asa = Wing(
        spans= [0.607431556265294, 0.5848349064082964],
        chords=[0.250693740050151, 0.250693740050151, 0.1253468700250755],
        offsets=[0, 0.0, 0.031336717506268874],
        twist_angles=[0, 0, 0],
        dihedral_angles=[0, 0, 0],
        airfoils= ['Marina03_alterado', 'Marina03_alterado', 'optfoilb2_alterado'],
        N_panels=12,
        distribution_type="cosine",
        sweep_check=False,
        surface_name='asa'
    )


asa.generate_mesh()

flight_condition = FlightCondition(
    V_inf=15,
    nu=1.3e-5,
    rho=1.086,
    aoa = [15, 17, 19],
    ground_effect_check=False
)

asa.setup_airfoil_data(flight_condition, airfoils_data_dict)

wingpool = WingPool(
    wing_list=[asa],
    flight_condition=flight_condition
)

simulation = Simulation(
    wing_pool=wingpool,
    damping_factor=0.5,
    max_iter=150,
    max_residual=1e-3,
    linear_check=False,
    show_logs=True,
    simulation_mode="latest_solution"
)

G_solution_list = simulation.run_simulation()

llt_coeficientes = PostProcessing(ref_point=[0, 0, 0])
CL_array = np.zeros(3)
CD_array = np.zeros(3)

for i, G in enumerate(G_solution_list):
    coefs = llt_coeficientes.get_global_coefficients(wingpool, G, aoa_index=i, S_ref=asa.total_area*2, c_ref=asa.MAC)
    CL_array[i] = coefs["CF"][2]
    CD_array[i] = coefs["CF"][0]

CLmax_dict = llt_coeficientes.get_CL_max_linear_new(wing_pool=wingpool, G_list=G_solution_list, S_ref=wingpool.S_ref)

CL_coeficientes = np.polyfit(flight_condition.aoa, CL_array, 1)
CD_coeficientes = np.polyfit(flight_condition.aoa, CD_array, 2)

asa_coeficientes = {
        "CLalfa": CL_coeficientes[0] * 180 / np.pi,
        "CL0": CL_coeficientes[1],
        "CDi_c1": CD_coeficientes[0] * ((180/np.pi)**2),
        "CDi_c2": CD_coeficientes[1] * (180/np.pi),
        "CDi_c3": CD_coeficientes[2],
        "Cmac": 0,
        "xCA": 0,
        "epsilon0": 180 / np.pi * 2 * CL_coeficientes[1] / (np.pi * asa.AR),
        "depsilondalpha": 180 / np.pi * 2 * CL_coeficientes[0] / (np.pi *  asa.AR)
    }

print(CL_array)

a=1