from ntt import *
import random
FINAL_DEGREE = 1  # Do not change it
from sympy import *
from Polynomial import *
from utilities import *

class BFV:
    N = None
    ntt_context = None
    PK = None
    SK = None

    def __init__(self, N, p, mu, sigma):
        self.N = N
        self.p = p
        self.mu = mu
        self.sigma = sigma

        self.q, k = generate_NTTPrime(2*N, 60)

        modulo_ntt = self.q
        zeta_ntt = PrimRoot_passing_from_generator(modulo_ntt, k)
        self.ntt_context = NTT_init(self.N, modulo_ntt, zeta_ntt, FINAL_DEGREE)

    def key_gen(self): # returns a secret key from R2
        self.SK = PolynomialGenerator.gen_from_uniform_distribution(2, self.q, self.N, self.ntt_context)
        PK2 = PolynomialGenerator.gen_from_uniform_distribution(self.q, self.q, self.N, self.ntt_context)
        e = PolynomialGenerator.gen_from_gaussian_distrribution(0, self.q, self.N, self.ntt_context, self.mu, self.sigma) # error with small coefficients

        PK1 = -(PK2*self.SK + e)
        self.PK = (PK1, PK2)

    def Relinearization_KeyGen(self, kk):
        a0 = PolynomialGenerator.gen_from_uniform_distribution(self.q*kk, self.q*kk, self.N, self.ntt_context)
        e = PolynomialGenerator.gen_from_gaussian_distrribution(0, self.q*kk, self.N, self.ntt_context, self.mu, self.sigma) # error with small coefficients
        
        Sk_1 = self.SK.to_R()
        Sk_2 = Sk_1*Sk_1
        RK1_0 = (a0.to_R()*Sk_1 + e.to_R())
        RK1_1 = (Sk_2*self.kk)
        RK1 = RK1_1 -RK1_0
        RK2 = a0.to_R_module(self.q*kk, self.ntt_context)
        RK = (RK1.to_R_module(self.q*kk, self.ntt_context), RK2)
        return RK

    def Relinearization_KeyGen_with_base_decomposition(self, base):
        l = int(floor(log(self.q,base)))
        SK_square = self.SK * self.SK

        RK = []
        # generating a0 and e for each component of the base decomposition
        for i in range(l+1):
            a0_i = PolynomialGenerator.gen_from_uniform_distribution(self.q, self.q, self.N, self.ntt_context)
            e_i = PolynomialGenerator.gen_from_gaussian_distrribution(0, self.q, self.N, self.ntt_context, self.mu, self.sigma) # error with small coefficients

            temp_sk_square = SK_square*(base**i)

            rk_i0 = -(a0_i*self.SK + e_i) + temp_sk_square
            rk_i1 = a0_i

            rk_i = (rk_i0, rk_i1)
            RK.append(rk_i)
        
        return RK

    def encrypt(self, M: Polynomial_ring):
        u = PolynomialGenerator.gen_from_uniform_distribution(2, self.q, self.N, self.ntt_context)
        e1 = PolynomialGenerator.gen_from_gaussian_distrribution(0, self.q, self.N, self.ntt_context, self.mu, self.sigma) # error with small coefficients
        e2 = PolynomialGenerator.gen_from_gaussian_distrribution(0, self.q, self.N, self.ntt_context, self.mu, self.sigma) # error with small coefficients
        PK1, PK2 = self.PK

        delta = int( self.q // self.p )
        a =  PK1*u + e1  + (M % self.q)*delta
        b = PK2*u + e2
        C = (a,b)
        return C

    def decrypt(self, C): # C is a couple of elements in Rq
        if len(C) == 2:
            a, b = C    
            delta = self.p/self.q
            return round((b*self.SK + a).to_R()*delta).to_R_module(self.p, self.ntt_context)
        elif len(C) == 3:
            C1, C2, C3 = C
            delta = self.p/self.q
            return round((C1 + C2*self.SK + C3*self.SK*self.SK).to_R()*delta).to_R_module(self.p, self.ntt_context)
         

    def homomorphic_addition(self,C1, C2): # C1,C2 is a couple of elements in Rq
        a1, b1 = C1
        a2, b2 = C2
        a = a1 + a2
        b = b1 + b2

        C = (a,b)
        return C

    def naive_homomorphic_multiplication(self,C1, C2):
        a1, b1 = C1
        a2, b2 = C2

        a1_tilde, a2_tilde, b1_tilde, b2_tilde = a1.to_R(), a2.to_R(), b1.to_R(), b2.to_R()
        delta = self.p / self.q

        C_1 = a1_tilde*a2_tilde
        C_2 = a1_tilde*b2_tilde + a2_tilde*b1_tilde
        C_3 = b1_tilde*b2_tilde
       
        C_1 = round(C_1*delta).to_R_module(self.q, self.ntt_context)
        C_2 = round(C_2*delta).to_R_module(self.q, self.ntt_context)
        C_3 = round(C_3*delta).to_R_module(self.q, self.ntt_context)
       
        C = (C_1, C_2, C_3)
        return C


    def homomorphic_multiplication(self,C1, C2, base_decomposition=True): # it uses Relinearization
        if base_decomposition:
            base = 256
            RK = self.Relinearization_KeyGen_with_base_decomposition(base)
            C_1, C_2, C_3 = self.naive_homomorphic_multiplication(C1, C2)

            l = int(floor(log(self.q,base)))
            C_3_decomposed = C_3.base_decomposition(base)

            C_3_0 = PolynomialGenerator.zero_polynomial(self.q, self.N, self.ntt_context)
            C_3_1 = PolynomialGenerator.zero_polynomial(self.q, self.N, self.ntt_context)
            for i in range(l+1):
                rk_i0, rk_i1 = RK[i]
                C_3_0 = C_3_0 + (rk_i0 * C_3_decomposed[i])
                C_3_1 = C_3_1 + (rk_i1 * C_3_decomposed[i])

            a = C_3_0 + C_1
            b = C_3_1 + C_2
            C = (a,b)
            return C
        else:
            kk = self.q**3 + 1
            RK1, RK2 = self.Relinearization_KeyGen(kk)
            C_1, C_2, C_3 = self.naive_homomorphic_multiplication(C1, C2)

            a = round((C_3.to_R()*RK1.to_R())/self.kk).to_R_module(self.q, self.ntt_context) + C_1
            b = round((C_3.to_R()*RK2.to_R())/self.kk).to_R_module(self.q, self.ntt_context) + C_2

            C = (a,b)
            return C
        

#########################################
#                 TEST                  #
#########################################

### bfv params
p = 32  # p should be a power of 2
N = 1024 # the degree of the polynomials will be N-1

# gaussian distribution
mu    = 0
sigma = 0.5 * 3.2

# bfv initialization
bfv = BFV(N, p, mu, sigma)
bfv.key_gen()

###### Creating two messages
m1, m2 = randint(0,2**(p-1)-1), randint(0,2**(p-1)-1)

M_1 = binary_encode_message(m1, N)
M_2 = binary_encode_message(m2, N)
M1 = Polynomial_ring(p, Polynomial(N, M_1), bfv.ntt_context) 
M2 = Polynomial_ring(p, Polynomial(N, M_2), bfv.ntt_context) 

print(f'M1 = {m1}, M2 = {m2}')

# Checking encryption and decryption
m1_after_decryption = decode_message(bfv.decrypt(bfv.encrypt(M1)).poly.coefficients, p)
assert m1 == m1_after_decryption


# Checking homomorphic addition
C1 = bfv.encrypt(M1)
C2 = bfv.encrypt(M2)
C3 = bfv.homomorphic_addition(C1, C2)
sum = m1 + m2
print(f'M1 + M2 = {sum}')
sum_after_homomorphic_addition = bfv.decrypt(C3)
assert sum == decode_message(sum_after_homomorphic_addition.poly.coefficients, p)

# Checking naive homomorphic multiplication (with SK^2)
multiplication = m1*m2
print(f'M1 * M2 = {multiplication}')
C4 = bfv.naive_homomorphic_multiplication(C1, C2)
multiplication_after_naive_homomorphic_mult = bfv.decrypt(C4)
assert multiplication == decode_message(multiplication_after_naive_homomorphic_mult.poly.coefficients, p)

# Checking homomorphic multiplication with RELINEARIZATION, without BASE DECOMPOSITION
# it does not work because the value C_3 * e increases the noise making the decryption infeasible
"""C5 = bfv.homomorphic_multiplication(C1, C2, base_decomposition=False)
multiplication_after_relinearization_homomorphic_mult = bfv.decrypt(C5)
assert multiplication == decode_message(multiplication_after_relinearization_homomorphic_mult.poly.coefficients, p)"""

# Checking homomorphic multiplication with RELINEARIZATION AND BASE DECOMPOSITION
C6 = bfv.homomorphic_multiplication(C1, C2, base_decomposition=True)
multiplication_after_relinearization_base_decomp_homomorphic_mult = bfv.decrypt(C6)
assert multiplication == decode_message(multiplication_after_relinearization_base_decomp_homomorphic_mult.poly.coefficients, p)

