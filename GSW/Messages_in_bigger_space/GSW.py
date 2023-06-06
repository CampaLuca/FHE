from MathObj.Matrices.matrix import *
from GSW.Messages_in_bigger_space.utilities import *
from MathObj.Number.Numbers import *
from sympy import log, sqrt, ceiling
from copy import deepcopy
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
        
        self.ZmodQ = Zmod(self.q)
        
        self.debug = debug
        if self.debug:
            print("Parameters:\n")
            print(f"k = {self.k}")
            print(f"n = {self.n}")
            print(f"q = {self.q}")

    def key_gen(self): # returns a secret key from R2
        
        ###### Generating the secret key
        self.s = Random_Matrix_Generator.from_uniform_distribution(self.q, self.n-1, 1, field=self.ZmodQ)
        self.sk = self.s.append_row([1])

        ###### generating the vector e (error)
        check = True
        while check:
            self.e = Random_Matrix_Generator.from_gaussian_distribution(2, self.mu, self.sigma, self.m, 1, field=self.ZmodQ)
            val = self.e.l_1_norm()
            if val < int(round(self.q/4)):
                check = False


        ###### generating the matrix A and the public key
        self.A = Random_Matrix_Generator.from_uniform_distribution(self.q, self.n-1, self.m, field=self.ZmodQ)
        self.pk = (-self.A).append_row( (self.s.transpose())*self.A + self.e.transpose() )

        ###### check the correctness of the generated key

        assert self.sk.transpose()*self.pk == self.e.transpose()

        # generating the Gadget Matrix (for base decomposition)
        self.G = generate_gadget_matrix(self.n, self.l)
        
        if self.debug:
            print("\nKeys:")
            print(f"Secret Key: {self.sk}")
            print(f"Public Key: {self.pk}")
            print(f"Error: {self.e}")

    """
        Input: a scalar in Z_
    """
    def encrypt(self, scalar):
        R = Random_Matrix_Generator.from_bernoulli_distribution(2, self.m, self.m, field=Zmod(2))
        C = self.pk * (R % self.q) + (scalar*self.G) % self.q
        return C

    def decrypt(self, C):
        numerator = (self.sk.transpose())*C
        denominator = (self.sk.transpose())*(self.G % self.q)

        assert numerator.get_dimension() == denominator.get_dimension()
        rows, columns = numerator.get_dimension()

        values = set()
        
        for c in range(columns):
            values.add( int(round(numerator[0,c].value / denominator[0,c].value, 0)) )
        
        minor_distance = float('inf')
        value = None

        for v in values:
            distance_matrix = (numerator - v*denominator) % self.q
            coefficients = distance_matrix.vector[0]
            coefficients2 = deepcopy(coefficients)
            for i in range(len(coefficients2)):
                coefficients2[i] = coefficients2[i].value
                
            for i in range(len(coefficients)):
                coefficients[i] = min(coefficients[i].value, self.q - coefficients2[i])

            distance = 0
            for a in coefficients:
                distance += a*a

            if distance < minor_distance:
                value = v
                minor_distance = distance
        
        return value


    """
        Given the encryption of two messages it returns the encryption of their sum in Z_q.
    """
    def homomorphic_addition(self, C1, C2):
        return C1 + C2

    """
        Given the encryption of two messages it returns the encryption of their product in Z_q.
    """
    def homomorphic_multiplication(self, C1, C2):
        C2_tilde = base_decomposed_matrix(C2, 2, self.n, self.l) % self.q # from n x m to m x m matrix
        return C1*C2_tilde
             

