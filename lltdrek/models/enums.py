from enum import Enum, EnumMeta

class CustomEnumMeta(EnumMeta):
    def __repr__(self):
        return ", ".join([f"{name}: {value}" for name, value in self.__members__.items()])

class DistributionTypes(Enum, metaclass=CustomEnumMeta):
    Linear = "linear"
    Cosine = "cosine"

class SimulationModes(Enum, metaclass=CustomEnumMeta):
    LinearFirst = "linear_first"
    LatestSolution = "latest_solution"
