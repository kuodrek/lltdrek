from lltdrek import Wing, FlightCondition, WingPool, Simulation, PostProcessing, load_folder

airfoils_data, _ = load_folder("sample_airfoils")
print(airfoils_data)

wing = Wing(
    spans=[1, 0.5],
    chords=[1, 0.5, 0.4],
    offsets=[0, 0, 0],
    twist_angles=[0, 0, 0],
    dihedral_angles=[0, 0, 0],
    airfoils=["NACA4424", "NACA4424", "NACA4412"],
    surface_name="wing",
    N_panels=4,
    distribution_type="cosine",
    sweep_check=False,
)
wing.generate_mesh()  # Generate wing simulation elements
flight_condition = FlightCondition(
    V_inf=60,
    nu=1.5e-5,
    rho=1.225,
    angles_of_attack=[0],
    h=0,
    ground_effect_check=False,
)
wing.setup_airfoil_data(flight_condition, airfoils_data)  # Assign airfoil data to wing

# Link wing that will be simulated with flight condition
wing_pool = WingPool(wing_list=[wing], flight_condition=flight_condition)

simulation = Simulation(
    damping_factor = 0.5,
    max_iter = 100,
    max_residual = 1e-4,
    linear_check = False,
    show_logs = False,
    simulation_mode = "linear_first" # Solves linear version of equations before doing non linear simulation to speed up process
)

# Run simulation

simulation_results = simulation.run(wing_pool)
import matplotlib.pyplot as plt

post_processing = PostProcessing()
coefficients = post_processing.get_coefficients(wing_pool, simulation_results)

CL_list = []
CM_list = []
CD_list = []

for coef in coefficients:
    CL_list.append(coef.global_coefficients.forces.CL)
    CM_list.append(coef.global_coefficients.moments.Cm)
    CD_list.append(coef.global_coefficients.forces.CD)


plt.scatter(flight_condition.angles_of_attack, CL_list)
plt.ylabel("CL")
plt.xlabel("Angle of attack")