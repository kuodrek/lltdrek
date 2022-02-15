'''
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
'''

class Wing():
    def __init__(self, b, c, offsets, tw_angles, dihedrals, airfoils, N_panels = 20, panels_distr = 'linear'):
        Wing.b = b
        Wing.c = c
        Wing.offsets = offsets
        Wing.twangles = tw_angles
        Wing.dihedrals = dihedrals
        Wing.airfoils = airfoils
        Wing.N_panels = N_panels
        Wing.panels_distr = panels_distr


class Flightconditions():
    def __init__(self, Vinf, nu, rho, h = 0):
            Flightconditions.Vinf = Vinf
            Flightconditions.nu = nu
            Flightconditions.rho = rho
            Flightconditions.h = h


class Simulation():
    def __init__(self, damp = 0.7, max_iter = 1000, convergence = 0.001):
        Simulation.damp = damp
        Simulation.max_iter = max_iter
        Simulation.convergence = convergence
