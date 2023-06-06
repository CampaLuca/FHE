import sys
sys.path.insert(0, '../')

from MathObj.Number.NTT import *

VERBOSE = False
FINAL_DEGREE = 1  # Do not change it

def test_NTT(ctxt, test_count):
    N = ctxt.N; q = ctxt.q
    f = [0]*N
    check = 0
    for j in range(test_count):
        for i in range(0, N):
            f[i]=random.randint(0,q-1)
        f_ = ntt(f, ctxt, FINAL_DEGREE) # compute NTT of f
        f__ = intt(f_, ctxt, FINAL_DEGREE) # compute inverse NTT of f_
        if poly_comp(f, f__)==False:
            check += 1 
    return check

def test_ringmult(ctxt, test_count):
    N = ctxt.N; q = ctxt.q
    f1 = [0]*N; f2 = [0]*N
    check = 0
    for j in range(test_count):
        for i in range(0, N):
            f1[i]=random.randint(0,q-1)
            f2[i]=random.randint(0,q-1)
        g1 = modulo_ringmul(f1, f2, q)
        g2 = ringmul_NTT(f1, f2, ctxt, FINAL_DEGREE)    
        if poly_comp(g1, g2) == False:
            check +=1
    return check

#random.seed(5)  # comment it out if you want randoms to change 

N = 128
q, k = generate_NTTPrime(2*N, 60)
zeta = PrimRoot_passing_from_generator(q, k)

input("Sayi gir")
  




ctxt = NTT_init(N, q, zeta, FINAL_DEGREE)

print("Testing NTT: Good if 0 --> ", test_NTT(ctxt, 2)) # works fine if it returns 0
print("Testing NTT based Ring Multiplication: Good if 0 --> ",test_ringmult(ctxt, 2)) # works fine if it returns 0 
