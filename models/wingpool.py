from typing import Dict, List
from models.wing import Wing
from models.flight_condition import FlightCondition
from utils import velocity
from dataclasses import dataclass, field
import numpy as np
import numpy.linalg as npla
import copy

'''
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
'''


@dataclass
class WingPool:
    wing_list: List[Wing]
    flight_condition: FlightCondition
    initial_G: Dict = None
    ind_velocities_dict: Dict = field(init=False)
    G_dict: Dict = field(init=False)
    total_velocity_dict: Dict = field(init=False)
    aoa_eff_dict: Dict = field(init=False)
    complete_wing_pool: List[Wing] = field(init=False)
    ind_velocities_list: List[Dict] = field(init=False)

    def __post_init__(self):
        self.ind_velocities_dict = {}
        self.G_dict = {}
        self.total_velocity_dict = {}
        self.aoa_eff_dict = {}
        self.complete_wing_pool = []
        self.ind_velocities_list = []

        self.build_complete_wing_pool()
        for _, wing in enumerate(self.complete_wing_pool):
            self.ind_velocities_dict[wing.surface_name] = {}
            self.total_velocity_dict[wing.surface_name] = []
            self.aoa_eff_dict[wing.surface_name] = {}
            if self.initial_G == None:
                self.G_dict[wing.surface_name] = [1 for _ in range(wing.N_panels)]
            else:
                for surface_name, G_list in self.initial_G.items():
                    self.G_dict[surface_name] = G_list
                    self.G_dict[surface_name+"_mirrored"] = G_list
        
        for v_inf_array in self.flight_condition.v_inf_array:
            ind_velocities_dict = self.calculate_induced_velocities(v_inf_array)
            self.ind_velocities_list.append(ind_velocities_dict)

    
    def build_complete_wing_pool(self):
        '''
        Method that builds a wing pool with mirrored wing objects.
        This is the list that will be used in calculations
        '''
        for wing in self.wing_list:
            mirrored_wing = copy.copy(wing)
            mirrored_wing.surface_name += "_mirrored"

            # Mirror y-coordinate of required values
            # TODO: do these operations using numpy logic
            for i in range(mirrored_wing.N_panels):
                mirrored_wing.u_a[i][1] *= -1
                mirrored_wing.u_n[i][1] *= -1
                mirrored_wing.collocation_points[i][1] *= -1
                mirrored_wing.vertice_points[i][1] *= -1
                mirrored_wing.cp_lengths[i][1] *= -1
                mirrored_wing.cp_dsl[i][1] *= -1
            mirrored_wing.vertice_points[-1][1] *= -1

            self.complete_wing_pool.append(wing)
            self.complete_wing_pool.append(mirrored_wing)

    def update_solution(self, G_solution) -> None:
        '''
        Method that splits the G_solution list, obtained by solving the 
        the main equations, into separate solution lists for each wing in the
        complete_wing_pool
        '''
        global_counter = 0
        for wing in self.complete_wing_pool:
            counter = 0
            G_list = []
            while counter < wing.N_panels:
                G_list.append(G_solution[global_counter])
                counter += 1
                global_counter += 1
            self.G_dict[wing.surface_name] = G_list
    
    def calculate_aoa_eff(self, total_velocity_dict: dict) -> dict:
        '''
        calcular o angulo de ataque efetivo de cada painel de cada asa
        a ideia é utilizar a lógica do numpy pra acelerar os cálculos
        precisa ser verificado se isso tá funcionando
        '''
        for wing in self.wing_list:
            self.aoa_eff_dict[wing.surface_name] = np.arctan(np.dot(total_velocity_dict[wing.surface_name], wing.u_n) / np.dot(total_velocity_dict[wing.surface_name], wing.u_a))
            self.aoa_eff_dict[wing.surface_name+"_mirrored"] = self.aoa_eff_dict[wing.surface_name]

    def calculate_induced_velocities(self, v_inf_array: np.ndarray) -> dict:
        # Essa função calcula tudo para um angulo de ataque e é chamada para cada aoa no __post_init__
        for wing_cp in self.complete_wing_pool:
            collocation_points = wing_cp.collocation_points
            cp_macs = wing_cp.cp_macs
            for wing_vp in self.complete_wing_pool:
                vertice_points = wing_vp.vertice_points
                v_ij_distr = velocity.get_induced_velocity_distribution(
                    collocation_points, cp_macs, vertice_points, v_inf_array
                )
                self.ind_velocities_dict[wing_cp.surface_name][wing_vp.surface_name] = v_ij_distr
        return self.ind_velocities_dict

    def calculate_total_velocity(self, v_inf_array: np.ndarray):
        '''
        WIP
        Method that calculates the sum of vij * G  + v_inf of all wings
        '''
        # É NECESSÁRIO ACHAR O INDICE DA DISTRIBUIÇAO DE VELOCIDADES CORRETO
        aoa_index = 0
        for idx, array in np.ndenumerate(self.flight_condition.v_inf_array):
            if array == v_inf_array:
                aoa_index = idx
        ind_velocities_dict = self.ind_velocities_list[aoa_index]

        total_velocity_panel = v_inf_array
        # Iterate through each wing
        for _, wing in enumerate(self.wing_list):
            total_velocity_distr = []
            G = self.G_dict[wing.surface_name]
            v_ij_distr = ind_velocities_dict[wing.surface_name]
            # Iterate through each panel of the target wing
            for _, panel in enumerate(v_ij_distr):
                # Iterate through each velocity induced by current wing
                for j, ind_velocity in enumerate(panel):
                    total_velocity_panel += ind_velocity * G[j]
                total_velocity_distr.append(total_velocity_panel)
            self.total_velocity_dict[wing.surface_name] = total_velocity_distr

        # validar se tá certo
        for wing_i in self.complete_wing_pool:
            for wing_j in self.complete_wing_pool:
                G = self.G_dict[wing_i.surface_name]
                v_ij_distr = ind_velocities_dict[wing_i.surface_name][wing_j.surface_name]
                total_velocity_distr = [] # velocidades de todos os paineis de uma asa
                for panel_i, v_ij_distr_panel in np.ndenumerate(v_ij_distr):
                    total_velocity_i = v_inf_array
                    # numpy logic: VALIDAR
                    # total_velocity_i += v_ij_distr[:]
                    for v_ij in v_ij_distr_panel:
                        total_velocity_i += v_ij
                    total_velocity_distr[panel_i] = total_velocity_i
            self.total_velocity_dict[wing_i.surface_name] = total_velocity_distr
                
