from scipy.spatial.distance import pdist, cdist, squareform

from ga_python.gapy import gago, bits2bytes
import numpy as np
from sklearn.metrics import pairwise_distances

gaoptions = {
    "PopulationSize": 300,
    "Generations": 100,
    "InitialPopulation": [],
    "MutationFcn": 0.05,
    "EliteCount": 2,
    "Verbose": False,
}

matriz_distancia = np.array([
    [0, 4, 8, 4, 4],
    [4, 0, 2, 9, 1],
    [8, 2, 0, 8, 5],
    [4, 9, 8, 0, 7],
    [4, 1, 5, 7, 0],
])

def fitness_function(bits, matrix=matriz_distancia):
    pontos = bits2bytes(bits)

    pontos_formatados = np.array(
        [pontos[0], pontos[1], pontos[2], pontos[3], pontos[4]],
        [pontos[5], pontos[6], pontos[7], pontos[8], pontos[9]],
    )

    distancias_estimadas = pdist(pontos_formatados, metric="euclidean")
    distancias_alvo = squareform(matrix)

    erro = distancias_estimadas - distancias_alvo

    soma_dos_erros_quadrados = np.sum(erro ** 2)

    max_abs_error = np.max(np.abs(erro))

    fitness = soma_dos_erros_quadrados + 5.0 * max_abs_error
    return fitness

results = gago(ffit=fitness_function, nbits=80, gaoptions=gaoptions)
print(results)