import numpy as np

"""
inputs iniciais
ASA
Vetor de envergaduras - b
vetor de cordas - c
vetor de offsets - offsets
vetor de twists - twangles
vetor de diedros - dihedral
vetor de perfis - airfoils
gap - para o caso de duas asas

CONDICOES DE VOO
Velocidade do escoamento - Vinf
Coeficiente de viscosidade dinamica do ar - nu
densidade do ar - rho
altura do chao - h
angulo de ataque - aoa
Numero de paineis - N
Distribuicao dos paineis - panel_distribution

SIMULACAO
Fator de amortecimento - damp
Fator de convergencia - R
Numero de iteracoes - max_iter
"""


class Wing:
    def __init__(
        self,
        spans,
        chords,
        offsets,
        twist_angles,
        dihedral_angles,
        airfoils,
        N_panels=20,
        distribution_type="linear",
    ):
        Wing.spans = spans
        Wing.chords = chords
        Wing.offsets = offsets
        Wing.twist_angles = [angle * np.pi / 180 for angle in twist_angles]
        Wing.dihedral_angles = [angle * np.pi / 180 for angle in dihedral_angles]
        Wing.airfoils = airfoils
        Wing.N_panels = N_panels
        Wing.distribution_type = distribution_type


class Flightconditions:
    def __init__(self, Vinf, nu, rho, h=0):
        Flightconditions.Vinf = Vinf
        Flightconditions.nu = nu
        Flightconditions.rho = rho
        Flightconditions.h = h


class Simulation:
    def __init__(self, damp=0.7, max_iter=1000, residual_value=0.001):
        Simulation.damp = damp
        Simulation.max_iter = max_iter
        Simulation.residual_value = residual_value
