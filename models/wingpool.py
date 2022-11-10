from typing import Dict, List
from models.wing import Wing
from models.flight_condition import FlightCondition
from utils import velocity
from dataclasses import dataclass, field
import numpy as np
import copy

'''
ind_velocities_dict = {
    "surface_name": {
        "surface_name": "v_ind_distr"
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
    initial_G: dict = None
    ind_velocities_dict: Dict = field(init=False)
    G_dict: Dict = field(init=False)
    total_velocity_dict: Dict = field(init=False)
    aoa_eff_dict: Dict = field(init=False)
    complete_wing_pool: List[Wing] = field(init=False)

    def __post_init__(self):
        self.ind_velocities_dict = {}
        self.G_dict = {}
        self.total_velocity_dict = {}
        self.aoa_eff_dict = {}
        self.complete_wing_pool = []

        self.build_complete_wing_pool()

        for _, wing in enumerate(self.complete_wing_pool):
            wing_name = wing.surface_name
            N_panels = wing.N_panels

            self.ind_velocities_dict[wing_name] = {}
            self.total_velocity_dict[wing_name] = []
            self.aoa_eff_dict[wing_name] = {}
            if self.initial_G == None:
                self.G_dict[wing_name] = []
                G = [1 for _ in range(N_panels)]
                self.update_solution(G=G,surface_name=wing_name)
            else:
                self.G_dict = self.initial_G

    def calculate_induced_velocities(self):
        v_inf = self.flight_condition.v_inf_array
        for wing_cp in self.complete_wing_pool:
            collocation_points = wing_cp.collocation_points
            cp_macs = wing_cp.cp_macs
            for wing_vp in self.complete_wing_pool:
                vertice_points = wing_vp.vertice_points
                v_ij_distr = velocity.get_induced_velocity_distribution(
                    collocation_points, cp_macs, vertice_points, v_inf
                )
                self.ind_velocities_dict[wing_cp.surface_name][wing_vp.surface_name] = v_ij_distr

        return self.ind_velocities_dict
    
    def build_complete_wing_pool(self):
        '''
        Method that builds a wing pool with mirrored wing objects.
        This list will be used in calculations
        '''
        for wing in self.wing_list:
            mirrored_wing = copy.copy(wing)
            mirrored_wing.surface_name += "_mirrored"

            # Mirror y-coordinate of required values
            for i in range(mirrored_wing.N_panels):
                mirrored_wing.u_a[1][i] *= -1
                mirrored_wing.u_n[1][i] *= -1
                mirrored_wing.collocation_points[1][i] *= -1
                mirrored_wing.vertice_points[1][i] *= -1
                mirrored_wing.cp_lengths[1][i] *= -1
                mirrored_wing.cp_dsl[1][i] *= -1
            mirrored_wing.vertice_points[1][-1] *= -1

            self.complete_wing_pool.append(wing)
            self.complete_wing_pool.append(mirrored_wing)

    def update_solution(self, G, surface_name: str):
        self.G_dict[surface_name] = G
    
    def calculate_aoa_eff(self):
        for wing in self.wing_list:
            total_velocity = self.total_velocity_dict[wing.surface_name]
            u_n_i = wing.u_n
            u_a_i = wing.u_a
            aoa_i = np.arctan(np.dot(total_velocity, u_n_i) / np.dot(total_velocity, u_a_i))
            a=1

    def get_total_dim_velocity(self, target_wing: Wing):
        '''
        Method that calculates the sum of vij * G  + v_inf of a target wing
        '''
        v_ij_dict = self.ind_velocities_dict[target_wing.surface_name]
        v_inf_array = self.flight_condition.v_inf_array
        v_total = v_inf_array

        # Iterate through each wing
        for _, wing in enumerate(self.wing_list):
            name = wing.surface_name
            G = self.G_dict[name]
            v_ij_distr = v_ij_dict[name]
            # Iterate through each panel of the target wing
            for _, panel in enumerate(v_ij_distr):
                # Iterate through each velocity induced by current wing
                for j, ind_velocity in enumerate(panel):
                    v_total += ind_velocity * G[j]

        return v_total

    def calculate_spanwise_distribution(self):
        '''
        TODO: solver1 equation
        '''
        pass

    def calculate_corrector_matrix(self):
        '''
        TODO: solver2 equation
        '''
        J_matrix = []
        for i, wing_i in enumerate(self.complete_wing_pool):
            row_i = []
            w_i = 1
            w_i_abs = 1
            v_n_i = 1
            v_a_i = 1
            C_l_alfa_i = 1
            for j, wing_j in enumerate(self.complete_wing_pool):
                # TODO: need to check when i = j, because a straight vortex induces no downwash on itself
                column_j = 0
                if wing_i.surface_name == wing_j.surface_name and i != j:
                    column_j += 2 * w_i_abs
            J_matrix.append(row_i)


