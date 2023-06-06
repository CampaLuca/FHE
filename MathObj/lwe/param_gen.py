from MathObj.Generators.prime_generators import *
from math import floor, log

def generate_parameters(k, mu, sigma):
    n = k
    q = generateSophieGermainPrime(k)
    l = floor(log(q, 2)) + 1
    m = n * l

    return k, n, m, l, q

def generate_parameters_for_gsw(k, mu, sigma, L):
    n = (k+110)/7.2
    q = 4
    l = floor(log(q,2)) + 1
    N = (n + 1) * l
    while (True):
        lower_bound = pow(N+1,L)
        lower_bound = lower_bound * 8 * int(sigma*6)
        if (q <= lower_bound):
            q = next_prime(q)
        else:
            break
        
        n = int(round(log(q/(ceiling(sigma)), 2)*(k+110)/(7.2)))
        l = floor(log(q,2)) + 1
        N = (n + 1) * l

    l = floor(log(q, 2)) + 1
    m = n * l

    return k, n, m, l, q
