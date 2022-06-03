import numpy as np
import math
import models

def generate_mesh(Wing: models.Wing):
    # Inicialização das variáveis
    N_panels = Wing.N_panels
    b_array = Wing.b
    b_abs = sum(b_array)
    c_array = Wing.c
    offsets_array = Wing.offsets
    dihedrals_array = Wing.dihedrals

    # Distribuição dos paineis por partição
    b_partitions = b_array / b_abs * N_panels
    b_partitions = [math.ceil(int(i)) for i in b_partitions]
    # Adicionar mais um painel por conta dos arredondamentos
    if sum(b_partitions) == N_panels - 1: b_partitions[0] += 1

    # Distribuição dos pontos de colocação e vértices
    N_partitions = len(b_array)
    collocation_points = np.array(np.zeros((N_panels, 3)))
    vertices = np.array(np.zeros((N_panels + 1, 3)))

    XYZ = np.array(np.zeros(N_partitions))
    XYZ_aux = [0, 0, 0]
    for i in range(N_partitions):
        XYZ[i][0] = 0.25*c_array[i+1] + offsets_array[i]
        XYZ[i][1] = XYZ_aux[1] + b_partitions[i] * np.cos(dihedrals_array * np.pi / 180)
        XYZ[i][2] = XYZ_aux[2] + np.tan(dihedrals_array * np.pi / 180) * b_partitions[i]
        XYZ_aux = XYZ[i][:]

    iter_auxX = 0
    iter_auxY = 0
    for i in range(N_partitions):
        if i == 1:
            X12 = [XYZ[i][0]- 0.25*c_array[i], XYZ[i][1], XYZ[i][2]]
            X1 = [0.25*c_array[i], 0, 0]
        else:
            X12 = [XYZ[i][0] - XYZ[i-1][0], XYZ[i][1] - XYZ[i-1][1], XYZ[i][2] - XYZ[i-1][2]]
            X1 = [XYZ[i-1][0], XYZ[i-1][1], XYZ[i-1][2]]
        
        if Wing.panels_distr == 'cosine':
            aspace1 = np.arange(0, np.pi, np.pi / b_partitions[i])
            aspace2 = np.arange(np.pi / (2 * b_partitions[i]), np.pi / b_partitions[i], np.pi - np.pi / (2 * b_partitions[i]))
            Lspace = (1 - np.cos(aspace1)) / 2
            Mspace = (1 - np.cos(aspace2)) / 2

            LspaceXZ = np.linspace(0, 1, b_partitions[i] + 1)
            MspaceXZ = np.zeros((1, b_partitions[i]))

            for j in range(aspace2):
                MspaceXZ[j] = (LspaceXZ[j] + LspaceXZ[j+1]) / 2

            X = np.zeros((b_partitions[i] + 1, 3))
            Y = np.zeros((b_partitions[i], 3))

            for k in range(aspace1):
                X[k][0] = X1[0] + LspaceXZ[k] * X12[0]
                X[k][1] = X1[1] + Lspace[k] * X12[1]
                X[k][2] = X1[2] + LspaceXZ[k] * X12[2]

            for k in range(aspace2):
                Y[k][0] = X1[0] + MspaceXZ[k] * X12[0]
                Y[k][1] = X1[1] + Mspace[k] * X12[1]
                Y[k][2] = X1[2] + MspaceXZ[k] * X12[2]

        elif Wing.panels_distr == 'linear':
            aspace1 = np.linspace(0, 1, b_partitions[i] + 1)
            aspace2 = np.zeros((1, b_partitions[i]))

            for j in range(aspace2):
                aspace2[j] = (aspace1[j] + aspace1[j + 1]) / 2
            
            Lspace = aspace1
            Mspace = aspace2

            X = np.zeros((b_partitions[i]+1, 3))
            Y = np.zeros((b_partitions[i], 3))

            for k in range(aspace1):
                X[k][0] = X1[0] + LspaceXZ[k] * X12[0]
                X[k][1] = X1[1] + Lspace[k] * X12[1]
                X[k][2] = X1[2] + LspaceXZ[k] * X12[2]

            for k in range(aspace2):
                Y[k][0] = X1[0] + MspaceXZ[k] * X12[0]
                Y[k][1] = X1[1] + Mspace[k] * X12[1]
                Y[k][2] = X1[2] + MspaceXZ[k] * X12[2]

        else:
            raise ValueError("Valor inválido em Wing.panels_distr: Esperava-se 'linear' ou 'cosine'")

        for j in range(X):
            vertices[iter_auxX][0] = X[j][0]
            vertices[iter_auxX][1] = X[j][1]
            vertices[iter_auxX][2] = X[j][2]
            iter_auxX += 1

        for j in range(Y):
            vertices[iter_auxY][0] = X[j][0]
            vertices[iter_auxY][1] = X[j][1]
            vertices[iter_auxY][2] = X[j][2]
            iter_auxY += 1






    print(b_partitions)
    print(N_partitions)

# Teste
b = np.array([3, 2, 1])
asa = models.Wing(b, 0, 0, 0, 0, 0)
generate_mesh(asa)