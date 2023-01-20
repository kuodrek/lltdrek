from typing import Dict, List
from models.wing import Wing
from models.flight_condition import FlightCondition
from utils import velocity
from dataclasses import dataclass, field
import numpy as np
import copy
from utils.timeit import timeit
from utils.lookup import get_airfoil_data


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


@dataclass(repr=False, eq=False, match_args=False)
class WingPool:
    wing_list: List[Wing]
    flight_condition: FlightCondition
    initial_G: Dict = None
    S_ref: float = None
    c_ref: float = None
    ind_velocities_dict: Dict = field(init=False)
    G_dict: Dict = field(init=False)
    complete_wing_pool: List[Wing] = field(init=False)
    total_panels: int = field(init=False)

    @timeit
    def __post_init__(self):
        self.G_dict = {}
        self.complete_wing_pool = []
        self.ind_velocities_list = []
        self.total_panels = 0
        self.S_ref
        self.c_ref 

        if self.S_ref is None:
            for _, wing in enumerate(self.complete_wing_pool):
                self.S_ref += wing.total_area * 2
        
        if self.c_ref is None:
            self.c_ref = self.wing_list[0].MAC
        
        if(len(set(self.wing_list)) == len(self.wing_list)):
            self.build_complete_wing_pool()
        else:
            raise Exception("As asas que compõem a wing pool precisam ter nomes únicos")
        
        for _, wing in enumerate(self.complete_wing_pool):
            self.total_panels += wing.N_panels
            if self.initial_G == None:
                self.G_dict[wing.surface_name] = [0.1 for _ in range(wing.N_panels)]
            else:
                for surface_name, G_list in self.initial_G.items():
                    self.G_dict[surface_name] = G_list
                    self.G_dict[surface_name+"_mirrored"] = G_list
        
        for v_inf_array in self.flight_condition.v_inf_list:
            ind_velocities_dict = self.calculate_induced_velocities(v_inf_array)
            self.ind_velocities_list.append(ind_velocities_dict)


    def build_complete_wing_pool(self):
        """
        Method that builds a wing pool with mirrored wing objects.
        This is the list that will be used in calculations
        """
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

            self.complete_wing_pool.append(wing)
            self.complete_wing_pool.append(mirrored_wing)


    def update_solution(self, G_solution) -> None:
        """
        Method that splits the G_solution list, obtained by solving the 
        the main equations, into separate solution lists for each wing in the
        complete_wing_pool
        """
        global_counter = 0
        for wing in self.complete_wing_pool:
            counter = 0
            G_array = np.zeros(wing.N_panels)
            while counter < wing.N_panels:
                G_array[counter] = G_solution[global_counter]
                counter += 1
                global_counter += 1
            self.G_dict[wing.surface_name] = G_array
        return self.G_dict
    
    @timeit
    def calculate_induced_velocities(self, v_inf_array: np.ndarray) -> dict:
        # Essa função calcula tudo para um angulo de ataque e é chamada para cada aoa no __post_init__
        ind_velocities_dict = {}
        for wing in self.complete_wing_pool:
            ind_velocities_dict[wing.surface_name] = {}
        
        for wing_cp in self.complete_wing_pool:
            collocation_points = wing_cp.collocation_points
            cp_macs = wing_cp.cp_macs
            for wing_vp in self.complete_wing_pool:
                vertice_points = wing_vp.vertice_points
                v_ij_distr = velocity.get_induced_velocity_distribution(
                    collocation_points, cp_macs, vertice_points, v_inf_array, wing_vp.surface_name, self.flight_condition.ground_effect_check, self.flight_condition.h
                )
                ind_velocities_dict[wing_cp.surface_name][wing_vp.surface_name] = v_ij_distr
        return ind_velocities_dict


    def calculate_total_velocity(self, v_inf_array: np.ndarray, G_dict: dict):
        """
        WIP
        Method that calculates the sum of vij * G  + v_inf of all wings
        - Dá pra acelerar esse método calculando somente a velocidade dos objetos originais e copiando para os objetos 
        espelhados
        """
        aoa_index = np.nan
        for idx, array in enumerate(self.flight_condition.v_inf_list):
            if np.array_equal(v_inf_array, array):
                aoa_index = idx
                break
        if aoa_index == np.nan:
            raise ValueError("v_inf_array inválido: o valor não se encontra na lista self.v_inf_list")
        
        ind_velocities_dict = self.ind_velocities_list[aoa_index]
        total_velocity_dict = {}
        for wing in self.complete_wing_pool:
            total_velocity_dict[wing.surface_name] = np.zeros((wing.N_panels,3))

        for wing_i in self.complete_wing_pool:
            for i, _ in enumerate(total_velocity_dict[wing_i.surface_name]):
                total_velocity_dict[wing_i.surface_name][i] += v_inf_array
                for wing_j in self.complete_wing_pool:
                    G = G_dict[wing_j.surface_name]
                    ind_velocities_distr = ind_velocities_dict[wing_i.surface_name][wing_j.surface_name][i]
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

    
    def get_Cl_distribution(self, aoa_eff_dict: dict, cl_alpha_check: bool) -> dict:
        """
        O método pode retornar tanto a distribuição de Cl quanto a distribuição de Clalpha
        do sistema de asas
        Não é utilizado atualmente
        """
        Cl_dict = {}
        for wing in self.wing_list:
            Cl_distr = np.zeros(wing.N_panels)
            aoa_eff_distr = aoa_eff_dict[wing.surface_name]
            for i, _ in enumerate(wing.collocation_points):
                Cl_distr[i] = get_airfoil_data(
                    wing.cp_airfoils[i],
                    wing.cp_reynolds[i],
                    aoa_eff_distr[i] * 180 / np.pi,
                    wing.airfoil_data,
                    cl_alpha_check = cl_alpha_check
                )
            Cl_dict[wing.surface_name] = Cl_distr
        return Cl_dict