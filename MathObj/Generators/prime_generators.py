import random

def next_prime(v):
    v = v + 1
    while not is_prime(v):
        v += 1
    
    return v

# This is not always true because there is a class of numbers whose property is to pass the fermat primality test
def is_prime(p):
    """ Returns whether p is probably prime """
    for null in range(16):
        a = random.randint(1,p-1)
        if pow(a,p-1,p) != 1:
            return False
    return True

# get Sophie Germain Prime (Safe prime of the form 2p + 1 with p prime)
def gen_prime(b):
    """ Returns a prime p with b bits """
    p = random.randint(2**(b-1), 2**b)
    while not is_prime(p):
        p = random.randint(2**(b-1), 2**b)
    return p

def generateSophieGermainPrime(k):
    """ Return a Sophie Germain prime p with k bits """
    p = gen_prime(k-1)
    sp = 2*p + 1
    while not is_prime(sp):
        p = gen_prime(k-1)
        sp = 2*p + 1
    return p

def generateSafePrime(k):
    """ Return a safe prime p with k bits """
    p = gen_prime(k-1)
    sp = 2*p + 1
    while not is_prime(sp):
        p = gen_prime(k-1)
        sp = 2*p + 1
    return sp