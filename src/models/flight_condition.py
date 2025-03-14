from typing import Union, Sequence
from src.models.types import AngleOfAttack
from src.models.exceptions import AlphaNotFoundException
import numpy as np

class FlightCondition:
    """Represents atmospheric conditions (air speed, density, etc) for a given
    wing pool

    :param V_inf: Free stream velocity
    :type V_inf: Union[float, int]
    :param nu: Air kinematic viscosity
    :type nu: Union[float, int]
    :param rho: Air density
    :type rho: Union[float, int]
    :param angles_of_attack: Sequence of angles of attack that will be used for simulation
    :type angles_of_attack: Sequence[AngleOfAttack]
    :param h: Height in reference to coordinates ```(x=0, z=0)``` in a wing pool system of surfaces.
    Used for ground effect calculations
    :type h: Union[float, int]
    :param ground_effect_check: If True, panel velocity calculation will take ground proximity into account
    :type ground_effect_check: bool = False
    """
    def __init__(
        self,
        V_inf: Union[float, int],
        nu: Union[float, int],
        rho: Union[float, int],
        angles_of_attack: Sequence[AngleOfAttack],
        h: Union[float, int],
        ground_effect_check: bool = False,
        angular_velocity: Sequence[float] = [0, 0, 0]
    ):
        self.V_inf = V_inf
        self.nu = nu
        self.rho = rho
        self.angles_of_attack = np.unique(angles_of_attack)
        self.h = h
        self.ground_effect_check = ground_effect_check
        self.v_inf_list = np.array([
            (
                np.cos(alpha * np.pi / 180),
                0,
                np.sin(alpha * np.pi / 180),
            )
            for alpha in self.angles_of_attack
        ])
        self.angular_velocity = np.array(angular_velocity)

    def get_alpha_index(self, alpha: AngleOfAttack) -> int:
        index = np.where(self.angles_of_attack == alpha)[0]
        if index.size == 0:
            raise AlphaNotFoundException(alpha, self.angles_of_attack)
        return index[0]

    def __repr__(self):
        return "FlightCondition({}={}, {}={}, {}={}, {}={}, {}={}, {}={})".format(
            "V_inf", self.V_inf, "nu", self.nu, "rho", self.rho, "angles_of_attack", self.angles_of_attack,
            "h", self.h, "ground_effect_check", self.ground_effect_check
        )
