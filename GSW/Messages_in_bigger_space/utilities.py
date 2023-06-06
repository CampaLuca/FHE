from MathObj.Matrices.matrix import *

# generate Gadget Matrix
def generate_gadget_matrix(n: int, l: int):
    M = Random_Matrix_Generator.zero_matrix(n, n*l)
    powers = [ 2**i for i in range(l)]

    for i in range(n):
        for j in range(l):
            M[i, i*l + j] = powers[j]

    return M

def base_decomposed_matrix(M: Matrix, base: int, n, l: int):
    rows, columns = M.get_dimension()
    output_matrix = Random_Matrix_Generator.zero_matrix(n*l, columns)
    
    if rows == n:
        for c in range(columns):
            for r in range(rows):
                value = M[r,c].value
                for b in range(l):
                    output_matrix[r*l + b, c] = (value & 0x1)
                    value = value >> 1
        
        return output_matrix
    else:
        raise Exception("Wrong input matrix dimension")
