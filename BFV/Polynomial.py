from random import randint,gauss
from sympy import *
from ntt import *


class Polynomial:
    def __init__(self, degree, coefficients):
        self.degree = degree
        self.coefficients = coefficients # from c_0 to c_(n-1)

class Polynomial_ring:
    def __init__(self, modulo, polynomial: Polynomial, ntt_ctxt):
        self.modulo = modulo
        self.poly = polynomial
        self.ntt_ctxt = ntt_ctxt

    def __add__(self, second_polynomial):
        # simple operation on ring without modulo
        if self.modulo is None and second_polynomial.modulo is None:
            new_poly = Polynomial(self.poly.degree, [(x+y) for x,y in zip(self.poly.coefficients, second_polynomial.poly.coefficients)])
            return Polynomial_ring(self.modulo, new_poly, self.ntt_ctxt)

        if self.modulo != second_polynomial.modulo:
            raise Exception("We are not in the same domain")

        new_poly = Polynomial(self.poly.degree, [(x+y) % self.modulo for x,y in zip(self.poly.coefficients, second_polynomial.poly.coefficients)])
        return Polynomial_ring(self.modulo, new_poly, self.ntt_ctxt)

    
    def __round__(self):
        coefficients = [int(round(c,0)) for c in self.poly.coefficients ]
        return Polynomial_ring(self.modulo, Polynomial(self.poly.degree, coefficients), self.ntt_ctxt)

    def __neg__(self):
        if self.modulo is None:
            coefficients = [(-c) for c in self.poly.coefficients ]
        else:
            coefficients = [(-c) % self.modulo for c in self.poly.coefficients ]

        return Polynomial_ring(self.modulo, Polynomial(self.poly.degree, coefficients), self.ntt_ctxt)
        
    
    def __sub__(self, second_polynomial):
        if self.modulo is None and second_polynomial.modulo is None:
            new_poly = Polynomial(self.poly.degree, [(x-y) for x,y in zip(self.poly.coefficients, second_polynomial.poly.coefficients)])
            return Polynomial_ring(self.modulo, new_poly, self.ntt_ctxt)

        if self.modulo != second_polynomial.modulo:
            raise Exception("We are not in the same domain")

        new_poly = Polynomial(self.poly.degree, [(x-y) % self.modulo for x,y in zip(self.poly.coefficients, second_polynomial.poly.coefficients)])
        return Polynomial_ring(self.modulo, new_poly, self.ntt_ctxt)

    def __mul__(self, value):
        if self.modulo is None and self.ntt_ctxt is None:
            if type(value) == float or type(value) == int:
                coefficients = [ (c*value) for c in self.poly.coefficients]
                return Polynomial_ring(self.modulo, Polynomial(self.poly.degree, coefficients), self.ntt_ctxt)
            elif type(value) == Polynomial_ring:
                coefficients = ringmul(self.poly.coefficients, value.poly.coefficients)
                return Polynomial_ring(self.modulo, Polynomial(self.poly.degree, coefficients), self.ntt_ctxt)
        elif self.modulo is not None and self.ntt_ctxt is not None:
            if type(value) == int:
                coefficients = [ (c*value) % self.modulo for c in self.poly.coefficients]
                return Polynomial_ring(self.modulo, Polynomial(self.poly.degree, coefficients), self.ntt_ctxt)
            elif type(value) == Polynomial_ring:
                if self.modulo != value.modulo:
                    raise Exception("We are not in the same domain")
                new_poly = Polynomial(self.poly.degree, ringmul_NTT(self.poly.coefficients, value.poly.coefficients, self.ntt_ctxt, 1))
                return Polynomial_ring(self.modulo, new_poly, self.ntt_ctxt)
        else:
            raise Exception("Unsupported operation")

    def __mod__(self, modulo):
        if self.modulo is not None and self.ntt_ctxt is not None:
            coefficients = [c % modulo for c in self.poly.coefficients ]
            return Polynomial_ring(modulo, Polynomial(self.poly.degree, coefficients), self.ntt_ctxt)
        else:
            raise Exception("Unsupported operation if polynomial coefficients are not integers")

    def __truediv__(self, value):
        if self.modulo is None and (type(value) == int or type(value) == float) and value != 0:
            coefficients = [ (c/value) for c in self.poly.coefficients]
            return Polynomial_ring(self.modulo, Polynomial(self.poly.degree, coefficients), self.ntt_ctxt)
        elif self.modulo is not None and type(value) == int and value != 0:
            try:
                inverse = modinv(value, self.modulo)
                coefficients = [ (c*inverse) % self.modulo for c in self.poly.coefficients]
                return Polynomial_ring(self.modulo, Polynomial(self.poly.degree, coefficients), self.ntt_ctxt)
            except:
                raise Exception(f"Modular Inverse of {value} does not exist under modulo {self.modulo}")
        else:
            raise Exception("None")

    def to_R_module(self, modulo, ntt_context):
        new_poly_ring = round(self)
        new_poly_ring.modulo = modulo
        new_poly_ring.ntt_ctxt = ntt_context
        return (new_poly_ring % modulo)

    def to_R(self):
        coefficients = [c for c in self.poly.coefficients ]
        return Polynomial_ring(None, Polynomial(self.poly.degree, coefficients), None)
    
    # decomposes each coefficient with respect to the supplied base
    # return a list of polynomials in ring
    def base_decomposition(self, base):
        l = int(floor(log(self.modulo,base)))
        decomposed_coefficients_polynomials_container = []

        coefficients = [c for c in self.poly.coefficients]

        for i in range(l+1):
            decomposed_coefficients_polynomials_container.append(PolynomialGenerator.zero_polynomial(self.modulo, self.poly.degree, self.ntt_ctxt))

        for i in range(l+1):
            for j in range(self.poly.degree):
                quotient = coefficients[j] // base
                remainder = coefficients[j] - quotient*base
                coefficients[j] = quotient
                decomposed_coefficients_polynomials_container[i].poly.coefficients[j] = remainder

        return decomposed_coefficients_polynomials_container



class PolynomialGenerator:
    def gen_from_uniform_distribution(B, modulo, N_values, ntt_ctxt):
        coefficients = [randint(-(B//2), B//2)% modulo for i in range(N_values)]

        new_poly = Polynomial(N_values, coefficients)
        return Polynomial_ring(modulo, new_poly, ntt_ctxt)
    
    def gen_from_gaussian_distrribution(B, modulo, N_values, ntt_ctxt, mu, sigma):
        coefficients = [int(gauss(mu,sigma))% modulo for i in range(N_values)]
        new_poly = Polynomial(N_values, coefficients)
        return Polynomial_ring(modulo, new_poly, ntt_ctxt)

    def zero_polynomial(modulo, N_values, ntt_ctxt):
        coefficients = [0 for i in range(N_values)]
        new_poly = Polynomial(N_values, coefficients)
        return Polynomial_ring(modulo, new_poly, ntt_ctxt)


    
    
    

