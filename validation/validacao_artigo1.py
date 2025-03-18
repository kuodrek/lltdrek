import numpy as np
from validacao_aoas import aoa_list_1
from lltdrek.models.wing import Wing
from lltdrek.models.flight_condition import FlightCondition
from lltdrek.models.wingpool import WingPool
from lltdrek.models.simulation import Simulation
from lltdrek.models.post_processing import PostProcessing
from lltdrek.utils import data


airfoils_data_dict, airfoils_dat_dict = data.load_folder('airfoils_test', -2, 6)

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
        distribution_type=DistributionTypes.Cosine,
        sweep_check=False,
        surface_name='asa'
    )


asa.generate_mesh()

flight_condition = FlightCondition(
    V_inf=45.48,
    nu=6.39e-6,
    rho=2.83,
    aoa = aoa_list_1,
    ground_effect_check=False
)

asa.setup_airfoil_data(flight_condition, airfoils_data_dict)

wingpool = WingPool(
    wing_list=[asa],
    flight_condition=flight_condition
)

simulation = Simulation(
    wing_pool=wingpool,
    damping_factor=1,
    max_iter=1000,
    max_residual=1e-3,
    linear_check=True,
    simulation_mode=SimulationModes.LinearFirst
)

G_solution_list = simulation.run_simulation()

llt_coeficientes = PostProcessing(ref_point=[0, 0, 0])

coefs_array = np.zeros((len(aoa_list_1), 3))

for i, _ in enumerate(G_solution_list):
    coefs = llt_coeficientes.get_global_coefficients(wingpool, G_solution_list[i], aoa_index=i, S_ref=asa.total_area*2, c_ref=asa.MAC)
    coefs_array[i,:] = np.array((coefs["CF"][0], coefs["CF"][2], coefs["CM"][1]))

np.savetxt("simulacao_linear.txt", coefs_array)
a=1