from typing import Dict, List, Sequence, Optional
import numpy as np
import copy

from src.models.wing import Wing
from src.models.flight_condition import FlightCondition
from src.utils import velocity
from src.utils.timeit import timeit
from src.utils.lookup import get_airfoil_data
from src.models.types import DVSArray, DVSMap
from src.models.exceptions import NonUniqueWingsException
from src.models.types import AngleOfAttack


"""
ind_velocities_list = [ind_velocities_dict_1, ind_velocities_dict_2, ...]

ind_velocities_dict = {
    "surface_name": {
        "surface_name": [[velocity_distribution_panel_1], [velocity_distribution_panel_2], ...]
    }
}

G_dict = {
    "surface_name": "G_array",
}

total_velocity_dict = {
   "surface_name": "total_velocity"
}

aoa_eff_dict = {
    "surface_name": "aoa_eff_list",
}
"""

class WingPool:
    def __init__(
        self,
        wing_list: List[Wing],
        flight_condition: FlightCondition,
        initial_G: Optional[Dict] = None,
        S_ref: Optional[float] = None,
        c_ref: Optional[float] = None,
        moment_ref: Sequence = [0, 0, 0]
    ):
        self.wing_list = wing_list
        self.pool: Sequence[Wing] = self._build_pool()
        self.flight_condition = flight_condition
        self.initial_G = initial_G
        self.S_ref = S_ref
        self.c_ref = c_ref
        self._moment_ref = np.array(moment_ref)
        self.system_moment_ref = self._build_system_moment_ref()

        self.G_dict = {}

        self.legacy_ind_velocities_list = []
        self.total_panels = 0

        if self.c_ref is None:
            self.c_ref = self.wing_list[0].MAC

        if self.S_ref is None:
            self.S_ref = sum([wing.total_area for wing in self.pool])

        for _, wing in enumerate(self.pool): # TODO: refactor this part of the code
            self.total_panels += wing.N_panels
            if self.initial_G == None:
                self.G_dict[wing.surface_name] = [0.1 for _ in range(wing.N_panels)]
            else:
                for surface_name, G_list in self.initial_G.items():
                    self.G_dict[surface_name] = G_list
                    self.G_dict[surface_name+"_mirrored"] = G_list
        
        self.system_freestream_velocities = self._build_system_freestream_velocities()
        self.system_induced_velocities = self._build_system_induced_velocities()

        for v_inf_array in self.flight_condition.v_inf_list: # TODO: switch this for loop to self.system_freestream_velocities
            ind_velocities_dict = self.legacy_calculate_induced_velocities(v_inf_array)
            self.legacy_ind_velocities_list.append(ind_velocities_dict)

    def _build_pool(self) -> List[Wing]:
        """Method that builds a wing pool with mirrored wing objects.
        This is the list that will be used in calculations
        """
        if not len(set(self.wing_list)) == len(self.wing_list):
            raise NonUniqueWingsException("All wings must have a unique name")

        wing_pool = []
        for wing in self.wing_list:
            mirrored_wing = copy.deepcopy(wing)
            mirrored_wing.surface_name += "_mirrored"

            # Mirror y-coordinate of required values
            mirrored_wing.u_a[:,1] = -1 * mirrored_wing.u_a[:,1]
            mirrored_wing.u_n[:,1] = -1 * mirrored_wing.u_n[:,1]
            mirrored_wing.u_s[:,0] = -1 * mirrored_wing.u_s[:,0]
            mirrored_wing.u_s[:,2] = -1 * mirrored_wing.u_s[:,2]
            mirrored_wing.collocation_points[:,1] = -1 * mirrored_wing.collocation_points[:,1]
            mirrored_wing.vertice_points[:,1] = -1 * mirrored_wing.vertice_points[:,1]
            mirrored_wing.cp_lengths[:,1] = -1 * mirrored_wing.cp_lengths[:,1]
            mirrored_wing.cp_dsl[:,0] = -1 * mirrored_wing.cp_dsl[:,0]
            mirrored_wing.cp_dsl[:,2] = -1 * mirrored_wing.cp_dsl[:,2]

            mirrored_wing.parent_wing = wing.surface_name

            wing_pool.append(wing)
            wing_pool.append(mirrored_wing)
        return wing_pool

    def map_solution(self, G: DVSArray) -> DVSMap:
        """Method that splits the DVS array, obtained by solving the 
        the main equations, into separate solution lists for each wing in
        wing pool
        """
        G_dict = {}
        global_counter = 0
        for wing in self.pool:
            counter = 0
            G_array = np.zeros(wing.N_panels)
            while counter < wing.N_panels:
                G_array[counter] = G[global_counter]
                counter += 1
                global_counter += 1
            G_dict[wing.surface_name] = G_array
        return G_dict

    def legacy_calculate_induced_velocities(self, v_inf_array: np.ndarray) -> dict:
        """**DEPRECATED**

        Pre calculates induced velocities for each panel and angle of attack
        """
        ind_velocities_dict = {}
        for wing in self.pool:
            ind_velocities_dict[wing.surface_name] = {}

        for wing_cp in self.pool:
            collocation_points = wing_cp.collocation_points
            cp_macs = wing_cp.cp_macs
            for wing_vp in self.pool:
                vertice_points = wing_vp.vertice_points
                v_ij_distr = velocity.get_induced_velocity_distribution(
                    collocation_points, cp_macs, vertice_points, np.tile(v_inf_array, (wing_cp.N_panels, 1)), wing_vp.surface_name, self.flight_condition.ground_effect_check, self.flight_condition.h
                )
                ind_velocities_dict[wing_cp.surface_name][wing_vp.surface_name] = v_ij_distr
        return ind_velocities_dict

    def calculate_total_velocity(self, alpha: AngleOfAttack, G_dict: dict):
        """
        Method that calculates the sum of vij * G  + v_inf of all wings
        - Dá pra acelerar esse método calculando somente a velocidade dos objetos originais e copiando para os objetos 
        espelhados
        - TODO: Test method with system_freestream_velocities
        """
        total_velocity_dict = {}
        for wing in self.pool:
            total_velocity_dict[wing.surface_name] = np.zeros((wing.N_panels,3))

        for wing_i in self.pool:
            wing_freestream_velocities = self.system_freestream_velocities[alpha][wing_i.surface_name]
            for i, _ in enumerate(total_velocity_dict[wing_i.surface_name]):
                total_velocity_dict[wing_i.surface_name][i] += wing_freestream_velocities[i]
                for wing_j in self.pool:
                    G = G_dict[wing_j.surface_name]
                    ind_velocities_distr = self.system_induced_velocities[alpha][wing_i.surface_name][wing_j.surface_name][i]
                    # ind_velocities_distr = ind_velocities_dict[wing_i.surface_name][wing_j.surface_name][i]
                    for j, v_ij in enumerate(ind_velocities_distr):
                        total_velocity_dict[wing_i.surface_name][i] += v_ij * G[j]
        return total_velocity_dict

    def calculate_aoa_eff(self, total_velocity_dict: dict) -> dict:
        """
        calcular o angulo de ataque efetivo de cada painel de cada asa
        """
        aoa_eff_dict = {}
        for wing in self.wing_list:
            total_velocity_distr = total_velocity_dict[wing.surface_name]
            aoa_eff_distr = np.zeros(wing.N_panels)
            for idx, velocity in enumerate(total_velocity_distr):
                aoa_eff_distr[idx] = np.arctan(np.dot(velocity, wing.u_n[idx]) / np.dot(velocity, wing.u_a[idx]))

            aoa_eff_dict[wing.surface_name] = aoa_eff_distr
            aoa_eff_dict[wing.surface_name+"_mirrored"] = aoa_eff_distr
        return aoa_eff_dict

    def _build_system_moment_ref(self) -> dict:
        system_moment_ref = {}
        for wing in self.pool:
            moment_ref_distribution = np.zeros((wing.N_panels, 3))
            moment_ref_distribution = wing.collocation_points - self._moment_ref
            system_moment_ref[wing.surface_name] = moment_ref_distribution
        return system_moment_ref

    def get_angular_velocity(self, angular_velocity) -> np.ndarray:
        system_angular_velocities = {}
        for wing, ref in self.system_moment_ref.items():
            system_angular_velocities[wing] = np.cross(angular_velocity, ref)
        return system_angular_velocities 
    
    def _build_system_freestream_velocities(self):
        # transformar 1 v_inf_array em um sistema de v_inf_arrays, um vetor tridimensional para cada painel (e isos para cada angulo)
        # retornar um dicionario onde a chave é o alfa e o valor é a distribuicao do sistema
        system_angular_velocity = self.get_angular_velocity(self.flight_condition.angular_velocity)
        system_freestream_velocities = {}
        for i, alpha in np.ndenumerate(self.flight_condition.angles_of_attack):
            freestream_array = self.flight_condition.v_inf_list[i]
            wing_freestream_velocities = {}
            for wing in self.pool:
                wing_freestream_velocities[wing.surface_name] = np.tile(freestream_array, (wing.N_panels, 1))  + system_angular_velocity[wing.surface_name]
            system_freestream_velocities[alpha] = wing_freestream_velocities

        return system_freestream_velocities

    def _build_system_induced_velocities(self):
        system_induced_velocities = {}
        for i, alpha in np.ndenumerate(self.flight_condition.angles_of_attack):
            freestream_velocities = self.system_freestream_velocities[alpha]
            wing_induced_velocities = {}
            for wing_cp in self.pool:
                wing_induced_velocities[wing_cp.surface_name] = {}

                wing_freestream_velocity = freestream_velocities[wing_cp.surface_name]
                collocation_points = wing_cp.collocation_points
                cp_macs = wing_cp.cp_macs
                for wing_vp in self.pool:
                    vertice_points = wing_vp.vertice_points
                    velocity_distribution = velocity.get_induced_velocity_distribution(
                        collocation_points, cp_macs, vertice_points, wing_freestream_velocity, wing_vp.surface_name, self.flight_condition.ground_effect_check, self.flight_condition.h
                    )
                    wing_induced_velocities[wing_cp.surface_name][wing_vp.surface_name] = velocity_distribution
            system_induced_velocities[alpha] = wing_induced_velocities

        return system_induced_velocities
