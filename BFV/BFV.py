from MathObj.Number.NTT import *
import random
FINAL_DEGREE = 1  # Do not change it
from sympy import *
from MathObj.Polynomials.Polynomial import *
from BFV.utilities import *

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
        
