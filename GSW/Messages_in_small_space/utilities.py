from MathObj.Matrices.matrix import *


# generate Gadget Matrix
def generate_gadget_matrix(n: int, l: int):
    M = Random_Matrix_Generator.zero_matrix(n, n*l)
    powers = [ 2**i for i in range(l)]

    for i in range(n):
        for j in range(l):
            M[i, i*l + j] = powers[j]

    return M

def base_decomposed_matrix(M: Matrix, base: int, l: int):
    rows, columns = M.get_dimension()
    output_matrix = Random_Matrix_Generator.zero_matrix(rows*l, columns)

    for c in range(columns):
        for r in range(rows):
            value = M[r,c].value
            for b in range(l):
                output_matrix[r*l + b, c] = (value & 0x1)
                value = value >> 1
    
    return output_matrix
    

def powers_of_base(M: Matrix, l: int, base: int):
    rows, columns = M.get_dimension()
    if columns > 1: 
        raise Exception("It is not a vector of dimension n x 1")

    N = Random_Matrix_Generator.zero_matrix(rows * l, 1)
    for i in range(rows):
        for exp in range(l):
            N[i*l + exp, 0] = M[i,0].value*2**exp

    return N


def flatten(M: Matrix, g: Matrix, l:int, base: int):
    return base_decomposed_matrix(g*M, base, l)

def LSB(v: int):
    return (v & 0x1)


    