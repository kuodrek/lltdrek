from typing import TypeAlias, Union, Dict
import numpy as np

AngleOfAttack: TypeAlias = Union[int, float]

# DVS Stands for Dimensionless Vortex Strenghts. These are the values we wish to find by solving the main equation
DVSArray: TypeAlias = np.ndarray # Unmapped array of DVS from Main Equation
DVSMap: TypeAlias = Dict[str, np.ndarray] # Map of DVS values, assigned to each wing in a wing pool
