import numpy as np

def calculate_theoretical_genotypes(p):
    q = 1 - p
    homo_p = round(p**2, 3)
    hetero = round(2 * p * q, 3)
    homo_q = round(q**2, 3)
    
    return homo_p, hetero, homo_q

def calculate_observed_genotypes(p, pop_size):
    population = (np.random.rand(pop_size, 2) < p).astype(int)
    
    homo_p, hetero, homo_q = 0, 0, 0

    pop_emoji = ""
    for i in range(pop_size):
        if np.sum(population[i,:]) == 2:
            pop_emoji += "🐦"
            homo_p += 1
        elif np.sum(population[i,:]) == 1:
            pop_emoji += "🦚"
            hetero += 1
        else:
            pop_emoji += "🐤"
            homo_q += 1
    
    homo_p /= pop_size
    hetero /= pop_size
    homo_q /= pop_size
    
    observed = homo_p, hetero, homo_q
    
    return pop_emoji, observed