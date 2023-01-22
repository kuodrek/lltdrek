from testes.validacao_aoas import aoa_list_1
from models.wing import Wing
from models.flight_condition import FlightCondition
from models.wingpool import WingPool
from models.simulation import Simulation
from models.post_processing import PostProcessing
from utils import data
import numpy as np

airfoils_data_dict, airfoils_dat_dict = data.load_folder('airfoils',aoa_polyfit_min=0,aoa_polyfit_max=10)

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
    # aoa = [aoa_list_1[0]],
    aoa = aoa_list_1,
    # aoa = [0, 10],
    ground_effect_check=False,
    h=0
)

asa.setup_airfoil_data(flight_condition, airfoils_data_dict)

wingpool = WingPool(
    wing_list=[asa],
    flight_condition=flight_condition
)

simulation = Simulation(
    wing_pool=wingpool,
    damping_factor=0.2,
    max_iter=50,
    max_residual=1e-3,
    linear_check=False,
    simulation_mode="linear_first",
    show_logs=True
)

G_solution_list = simulation.run_simulation()

llt_coeficientes = PostProcessing(ref_point=[0, 0, 0])

output_array = np.zeros((len(G_solution_list),3))

for idx, G in enumerate(G_solution_list):
    aoa = wingpool.flight_condition.aoa[idx]
    coefs = llt_coeficientes.get_global_coefficients(wingpool, G_solution_list[idx], aoa_index=idx, S_ref=wingpool.S_ref, c_ref=wingpool.c_ref)
    CF = coefs["CF"]
    CM = coefs["CM"]
    output_array[idx,:] = np.array((aoa, CF[2], CM[1]))
CLmax = llt_coeficientes.get_CL_max_linear(wingpool, G_solution_list, None, S_ref=asa.total_area*2)

# ac_data = llt_coeficientes.get_aerodynamic_center(wing_pool=wingpool, G_list=G_solution_list, aoa_1=0, aoa_2=10, S_ref=wingpool.S_ref, c_ref=wingpool.c_ref)

np.savetxt("output_final.txt", output_array)

a=1