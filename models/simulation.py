from dataclasses import dataclass


@dataclass
class Simulation:
    damp: float = 0.7
    max_iter: int = 1000
    residual_value: float = 0.001
    linear_check: bool = False

