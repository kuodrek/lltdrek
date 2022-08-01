class Simulation:
    """
    Fator de amortecimento - damp
    Fator de convergencia - R
    Numero de iteracoes - max_iter
    """
    def __init__(self, damp=0.7, max_iter=1000, residual_value=0.001):
        Simulation.damp = damp
        Simulation.max_iter = max_iter
        Simulation.residual_value = residual_value
