from sympy import log, sqrt, ceiling
from copy import deepcopy
from GSW.Messages_in_small_space import utilities
from MathObj.Matrices.matrix import *
from MathObj.Number.Numbers import *
from MathObj.Generators.prime_generators import *
from MathObj.lwe.param_gen import *


class GSW:
    N = None
    ntt_context = None
    PK = None
    SK = None

    def __init__(self, k, mu, sigma, debug=False): # k is the number of bits 
        
        self.n = k

        # parameters for gaussian distribution
        self.mu = mu
        self.sigma = sigma
        
        self.k, self.n, self.m, self.l, self.q = generate_parameters(k, mu, sigma)
        self.N = (self.n+1)*self.l  
        self.ZmodQ = Zmod(self.q)
        
        self.debug = debug
        if self.debug:
            print("Parameters:\n")
            print(f"k = {self.k}")
            print(f"n = {self.n}")
            print(f"q = {self.q}")
            

    def key_gen(self): 
        
        ###### Generating the secret key
        self.t = Random_Matrix_Generator.from_uniform_distribution(self.q, self.n, 1, field=self.ZmodQ)
        self.sk = (-self.t).prepend_row([1])
        self.v = utilities.powers_of_base(self.sk, self.l, 2)

       
        ###### Generating the error e
        check = True
        while check:
            self.e = Random_Matrix_Generator.from_gaussian_distribution(2, self.mu, self.sigma, self.m, 1, field=self.ZmodQ)
            val = self.e.l_1_norm()
            if val < int(floor(self.q/12)):
                check = False


        ###### generating the matrix A and the public key
        self.A = Random_Matrix_Generator.from_uniform_distribution(self.q, self.n, self.m, field=self.ZmodQ)

        self.pk = self.A.prepend_row( (self.t.transpose())*self.A + self.e.transpose() )

        # check the correctness of the generated key
        assert self.sk.transpose()*self.pk == self.e.transpose()

        # generating the Gadget Matrix (for base decomposition)
        self.G = utilities.generate_gadget_matrix(self.n+1, self.l) % self.q

        if self.debug:
            print("\nKeys:")
            print(f"Secret Key: {self.sk}")
            print(f"Public Key: {self.pk}")
            print(f"Error: {self.e}")


    """
        Input: a scalar in {0,1}
    """
    def encrypt(self, scalar):
        R = Random_Matrix_Generator.from_bernoulli_distribution(2, self.m, self.N, field=Zmod(2))
        I_N = Random_Matrix_Generator.identity_matrix(self.N, field=self.ZmodQ)
        
        part1 = utilities.base_decomposed_matrix(self.pk * (R % self.q), 2, self.l) % self.q
        part2 = (scalar*I_N) % self.q

        C = utilities.flatten((part1+part2), self.G, self.l, 2) % self.q
        return C

    def decrypt(self, C):
        index = -1

        for i in range(self.l):
            if self.v[i,0].value > self.q/4 and self.v[i,0].value <= int(floor(self.q/2)):
                index = i
                break

        if index == -1:
            raise Exception("We cannot find a vector which fits the property")
            
        C_i = C[None,index] % self.q
        C_out = C_i.element_wise_product(self.v % self.q)

        x = 0
        for r in range(C_out.get_dimension()[0]):
            x += C_out[r,0].value

        x = x % self.q
        v_i = self.v[index,0].value
        result = int(round(x/v_i))

        return result
        

    """
        Given the encryption of messages in {0,1} it returns the encryption of their sum in Z.
    """
    def homomorphic_addition(self, C1, C2):
        return utilities.flatten(C1 + C2 , self.G, self.l, 2) % self.q

    """
        Given the encryption of messages in {0,1} it returns the encryption of their product
    """
    def homomorphic_multiplication(self, C1, C2):
        return utilities.flatten(C1*C2, self.G, self.l, 2) % self.q

    """
        Given the encryption of messages in {0,1} it returns the encryption of the NAND value
    """
    def NAND(self, C1, C2):
        I_N = Random_Matrix_Generator.identity_matrix(self.N, self.ZmodQ)
        return utilities.flatten(I_N - C1*C2, self.G, self.l, 2) % self.q
             
