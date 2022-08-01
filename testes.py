from models import Wing

asa = Wing(
    spans=[3, 2],
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
print(asa.collocation_points)