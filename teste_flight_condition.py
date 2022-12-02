from models.flight_condition import FlightCondition


flight_condition = FlightCondition(
    V_inf=15,
    nu=1.5e-5,
    rho=1.225,
    aoa = [1, 2, 3],
    ground_effect_check=False
)

a=1