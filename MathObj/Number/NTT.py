##### Credits to: Erkay Savas from Sabanci University


import math
import time
import random
import sys
import sympy

DEBUG = False
DOUBLEROU = False

#########################################################
# This block generates multiplicative inverse
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('Modular inverse does not exist')
    else:
        return x % m
#########################################################

def generate_NTTPrime(N, bit_length, ):
    while True:
        p = random.randint(0,2**bit_length)
        if DOUBLEROU:
            if sympy.isprime(p) == True and (p-1)%(2*N) == 0:
                break
        else:
            if sympy.isprime(p) == True and (p-1)%N == 0 and (p-1)%(2*N) != 0:
                break
    
    if DOUBLEROU:
        for k in range(bit_length-1, 0, -1):
            if (p-1) % pow(2,k) == 0:
                return p, (p-1)//pow(2,k)
        return None
    else:
        return p, ((p-1)//N)


# factorization
def get_factors(n):
    return list(sympy.factorint(n))


# Primitive root of unity hesaplanmasinda kullanilabilir.
def PrimRoot_old(q, m):
    while True:
        zeta = random.randint(0,q-1)
        chck = True
        for i in range(1, m):
            if pow(zeta, i, q) == 1:
                chck = False
                break
        if pow(zeta, m, q) == 1 and chck:
            return zeta
    return 0

def PrimRoot(q, m):
    while True:
        zeta = random.randint(0,q-1)
        if pow(zeta, m, q) == 1:
            if pow(zeta, m//2, q) == q-1:
                return zeta
    return 0

def PrimRoot_passing_from_generator(q, k):
    # if DOUBLEROT is True k = 2^k'
    # else k = integer
    if DOUBLEROU:
        factor = 2
        generator = None
        for g in range(1,q):
            check = (pow(g, (q-1)//2, q) != 1)            
            if check:
                generator = g
                break
        
        if generator is None:
            return None
        
        return pow(generator, (q-1)//k,q)
    else:
        generator = None
        factors = get_factors(q-1)
        for g in range(1,q):
            check = True
            for factor in factors:
                if pow(g, (q-1)//factor, q) == 1:
                    check = False
                    break
            
            if check:
                generator = g
                break
        
        if generator is None:
            return None
        #return generator
        return pow(generator,k,q) # w = g^k


def PrimRoot_directly(q, k):
    # if DOUBLEROT is True k = 2^k'
    # else k = integer
    n = ((q - 1) // k)
    list_of_factors = get_factors(n)

    for a in range(1,q):
        if pow(a,n,q) == 1:
            check = True
            for factor in list_of_factors:
                if pow(a,(n//factor),q) == 1:
                    check = False
                    break
            if check:
                return a
    
    return None
                    



def br(i, n):
    return int(format(i, '0%db' % n)[::-1], 2)

def gen_powers(N, q, zeta, findeg):
    if findeg == 1 or findeg == 2:
        if findeg == 1: N_ = N
        else: N_ = N//2
        powers = [0]*N_    
        powers[0] = 1
        powers[1] = N_//2
        i = 1
        while 2**i < N_:
            for j in range(2**i, 2**(i+1), 2):
                powers[j] = powers[j//2]//2
                powers[j+1] = (powers[j//2]+N_)//2
            i = i + 1
    elif findeg == 3:
        N_ = N//2
        powers = [0]*((N//3))    
        powers[0] = 0
        powers[1] = N//6
        powers[2] = powers[1]//2
        powers[3] = (5*powers[1])//2
        i = 2
        while 2**i < N//3:
            for j in range(2**i, 2**(i+1)-2**(i-1), 2):
                powers[j] = powers[j//2]//2
                powers[j+1] = (powers[j//2]+N_)//2
            for j in range(2**(i+1)-2**(i-1), 2**(i+1), 2):
                powers[j] = powers[j//2]//2
                powers[j+1] = (powers[j//2]+N_)//2    
            i = i + 1
    return powers

def gen_twiddles(N2, q, zeta, powers):
    twiddle_cnt = len(powers)
    twiddles = [0]*twiddle_cnt
    inv_twiddles = [0]*twiddle_cnt 
    tmp = [0]*(N2)
    for i in range(N2):
        tmp[i] = pow(zeta, i, q)
    for i in range(twiddle_cnt):
        twiddles[i] = tmp[powers[i]]
        inv_twiddles[i] = -(tmp[powers[twiddle_cnt-1-i]])%q    
    return twiddles, inv_twiddles            

def poly_comp(f1, f2):
    result = True
    for i in range(0, len(f1)):
        if f1[i]!=f2[i]:
            result = False
    return result        

class NTT_ctxt:
    N = 0
    q = 0
    zeta = 0
    zetas = []
    zetas_inv = []
    post_proc = 0
    NTTinv = 0
    prec = 0

    def __init__(self):
        N = 0
    

    
def NTT_init(N, q, zeta, findeg):
    ntt_ctxt = NTT_ctxt()
    ntt_ctxt.N = N
    ntt_ctxt.q = q
    ntt_ctxt.zeta = zeta
    if findeg == 1:
        N2 = N
    elif findeg == 2:
        N2 = N//2
    elif findeg == 3:
        N2 = N//2
    ntt_ctxt.prec = int(math.log(N2, 2))
    powers = gen_powers(N, q, zeta, findeg)
    ntt_ctxt.zetas, ntt_ctxt.zetas_inv = gen_twiddles(N2, q, zeta, powers)
    #for i in range(0, N2):
        #NTT_ctxt.zetas.append(pow(zeta, br(i, NTT_ctxt.prec), q))
        #NTT_ctxt.zetas_inv.append((-pow(zeta, N2-1-br(i, NTT_ctxt.prec), q))%q)
    ntt_ctxt.post_proc = modinv((2*ntt_ctxt.zetas[1]-1)%q, q) # (2*zeta^(N//2)-1)^(-1)
    if findeg == 3:
        ntt_ctxt.NTTinv = modinv(N//3, q)
    else:
        ntt_ctxt.NTTinv = modinv(N2, q)
    return ntt_ctxt

def ntt(ff, ntt_ctxt, findeg):
    f = [0]*len(ff)
    for i in range(0, len(ff)):
        f[i] = ff[i] 
    N = ntt_ctxt.N
    q = ntt_ctxt.q
    mult_cnt = 0
    if findeg == 3:
        length = N//4
        for i in range(N//2):
            omega = ntt_ctxt.zetas[1] 
            t = (omega*f[i+N//2])%q
            f[i+N//2] = (f[i]+f[i+N//2]-t)%q
            f[i] = (f[i]+t)%q
            mult_cnt += 1
        k = 2
    else:
        length = N//2
        k = 1
    j = 0    
    while length >= findeg:
        start = 0
        while start < N:
            if DEBUG: input("push any button")
            omega = ntt_ctxt.zetas[k]
            if DEBUG: print("zeta power: ", br(k, NTT_ctxt.prec))
            k += 1
            for j in range(start, start+length):
                if DEBUG: print("i, j: ", j, j+length)
                t = (omega*f[j+length])%q
                f[j+length] = (f[j] - t)%q
                f[j] = (f[j] + t)%q
                mult_cnt += 1
            start = (j+1) + length
        length = length//2
    #print("Mult count (NTT): ", mult_cnt)    
    return f    
        
def intt(ff, ntt_ctxt, findeg):
    f = [0]*len(ff)
    for i in range(0, len(ff)):
        f[i] = ff[i] 
    N = ntt_ctxt.N
    q = ntt_ctxt.q
    mult_cnt = 0
    j = 0
    k = 0
    length = findeg
    if findeg == 3:
        N_ = N//4
    else:
        N_ = N//2
    while length <= N_:
        start = 0
        while start < N:
            omega = ntt_ctxt.zetas_inv[k]
            k += 1
            for j in range(start, start+length):
                t = f[j]
                f[j] = (t+f[j + length])%q
                f[j + length] = (t - f[j + length])%q
                f[j + length] = (omega*f[j + length])%q
                mult_cnt += 1
            start = (j+1) + length
        length = (length << 1)

    if findeg == 3:
        for i in range(0, N//2):
            t = ((f[i]-f[i+N//2])*ntt_ctxt.post_proc)%q   # t = f[i+N//2]
            f[i] = (f[i]+f[i+N//2]-t)%q
            f[i+N//2] = (2*t)%q
    for j in range(0, N):
        f[j] = (f[j]*ntt_ctxt.NTTinv)%q
        mult_cnt += 1
    #print("Mult count (INTT): ", mult_cnt)        
    return f

#reference polynomial multiplication function
def modulo_ringmul(f1, f2, q):
    C = [0]*(2*len(f1))
    D = [0]*(len(f1))
    for index1,elem1 in enumerate(f1):
        for index2,elem2 in enumerate(f2):
            C[index1+index2] = (C[index1+index2] + elem1*elem2)%q

    for i in range (len(f1)):
        D[i] = (C[i] - C[i+len(f1)])%q

    return D

def ringmul(f1, f2):
    C = [0]*(2*len(f1))
    D = [0]*(len(f1))
    for index1,elem1 in enumerate(f1):
        for index2,elem2 in enumerate(f2):
            C[index1+index2] = (C[index1+index2] + elem1*elem2)

    for i in range (len(f1)):
        D[i] = (C[i] - C[i+len(f1)])

    return D

def ringmul_NTT(f1, f2, ntt_ctxt, degree):
    N = ntt_ctxt.N
    q = ntt_ctxt.q
    f1_ = ntt(f1, ntt_ctxt, degree)
    f2_ = ntt(f2, ntt_ctxt, degree)
    z_tilde = []
    for i in range(len(f1)):
        z_tilde.append( (f1_[i]*f2_[i] % q) )
    z = intt(z_tilde, ntt_ctxt, degree)
    return z
    
    
    
