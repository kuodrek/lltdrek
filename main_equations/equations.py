import numpy as np
from models.wing import Wing
from models.simulation import Simulation
from models.flight_conditions import FlightConditions

def main_equation(Wing: Wing, Simulation: Simulation, FlightConditions: FlightConditions):
    residual = 0
    return residual

def corrector_equation(Wing: Wing, Simulation: Simulation, FlightConditions: FlightConditions):
    delta_G = 0
    return delta_G
    