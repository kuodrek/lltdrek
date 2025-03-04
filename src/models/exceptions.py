# Flight condition exceptions
class AlphaNotFoundException(IndexError):
    def __init__(self, alpha, angles_of_attack):
        super().__init__(
            f"Angle of attack {alpha} not found in angles_of_attack array: {angles_of_attack}"
        )


# Wing Pool exceptions
class NonUniqueWingsException(Exception): ...
