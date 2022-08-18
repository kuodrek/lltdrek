import numpy as np
import timeit
from models.wing import Wing

N_panels = 4
mat_teste = np.zeros([N_panels, N_panels, 3])
G = np.zeros(N_panels)

for i in range(N_panels):
    G[i] = i
    for j in range(N_panels):
        v_ij = np.array([i+1, j+1, i+j+1])
        mat_teste[i,j,:] = v_ij

for i in range(N_panels):
    final = mat_teste[i,:,:] * G[i]

print(G)
print(mat_teste)
print(final)