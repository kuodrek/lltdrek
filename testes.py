import numpy as np
import timeit
from models.wing import Wing

N_iter = 100
start_time = timeit.default_timer()
for i in range(N_iter):
    asa = Wing(
        spans=[3, (1+i/N_iter)],
        chords=[1, 0.8, 0.4],
        offsets=[0, 0, 0.5],
        twist_angles=[0, 0, 0],
        dihedral_angles=[10, 15],
        airfoils=['optfoilb2', 'optfoilb2', 'optfoilb2'],
        N_panels=12,
        distribution_type="cosine",
        sweep_check=False
    )
    asa.generate_mesh()
    # print(f"Iteracao {i}")
print(f"Tempo total: {timeit.default_timer()-start_time}")
# np.show_config()