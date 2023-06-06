from copy import deepcopy
import random
from MathObj.Number.Numbers import *
from MathObj.Generators.discretegauss import sample_dgauss, sample_bernoulli_exp
from fractions import Fraction

class Matrix:
    def __init__(self, values: list, field, rows=None, columns=None): # values must be a well formatted matrix [ [ [], [], [], [] ] ]
        if len(values) == 0:
            raise Exception("Empty list")

        if rows is None or columns is None:
            self.field = field
            
            rows = len(values)
            columns = len(values[0])

            for r in range(1,rows):
                if len(values[r]) != len(values[r-1]):
                    raise Exception("The rows do not have the same length")

            
            for r in range(rows):
                for c in range(columns):
                    values[r][c] = field(values[r][c])
            self.vector = values
        else:
            if len(values) != rows*columns:
                raise Exception("You supplied less values then required")

            self.field = field
            self.vector = []
            
            for r in range(rows):
                temp_list = []
                for c in range(columns):
                    temp_list.append(field(values[r*columns + c]) )
                    

                self.vector.append(temp_list)

    
    def __str__(self):
        dim = self.get_dimension()
        result = ""
        for r in range(dim[0]):
            for c in range(dim[1]):
                result = result + str((self.vector[r][c]).value) + " "
            result += "\n"
        return result

    def get_dimension(self):
        if type(self.vector[0]) != list:
            return ( len(self.vector), 1 )
        return ( len(self.vector), len(self.vector[0]) ) 
    
    def transpose(self): # modifies the object
        rows, columns = self.get_dimension()
        temp_vector = []
        for c in range(columns):
            t = []
            for r in range(rows):
                t.append(0)
            temp_vector.append(t)

        for r in range(rows):
            for c in range(columns):
                temp_vector[c][r] = self.vector[r][c]
        
        return Matrix(temp_vector, self.field)
    

    def __add__(self, second_matrix):
        if self.field == second_matrix.field:
            self_dimensions = self.get_dimension()
            second_matrix_dimension = second_matrix.get_dimension()
            return_vector = deepcopy(self)
            if self_dimensions == second_matrix_dimension:
                for row in range(self_dimensions[0]):
                    for col in range(self_dimensions[1]):
                        return_vector.vector[row][col] = self.vector[row][col] + second_matrix.vector[row][col]

                return return_vector

            else:
                raise Exception(f"Cannot sum matrix of dimension {self_dimensions} with matrix of dimension {second_matrix_dimension}")
        else:
            raise Exception("They are not in the same domain")

    def __sub__(self, second_matrix):
        if self.field == second_matrix.field:
            self_dimensions = self.get_dimension()
            second_matrix_dimension = second_matrix.get_dimension()
            return_vector = deepcopy(self)
            if self_dimensions == second_matrix_dimension:
                for row in range(self_dimensions[0]):
                    for col in range(self_dimensions[1]):
                        return_vector.vector[row][col] = self.vector[row][col] - second_matrix.vector[row][col]

                return return_vector

            else:
                raise Exception(f"Cannot sum matrix of dimension {self_dimensions} with matrix of dimension {second_matrix_dimension}")
        else:
            raise Exception("They are not in the same domain")

    

    def __truediv__(self, scalar):
        if type(scalar) == int or type(scalar) == float:
            self_dimensions = self.get_dimension()
            temp_vector = []
            for r in range(self_dimensions[0]):
                t = []
                for c in range(self_dimensions[1]):
                    t.append(0)
                temp_vector.append(t)
            for r in range(self_dimensions[0]):
                for c in range(self_dimensions[1]):
                    temp_vector[r][c] = self.vector[r][c] / scalar

            return Matrix(temp_vector, field=self.field)
        else:
            raise Exception("Unsuppoted operation")

    def __neg__(self):
        self_dimensions = self.get_dimension()
        temp_vector = []
        for r in range(self_dimensions[0]):
            t = []
            for c in range(self_dimensions[1]):
                t.append(0)
            temp_vector.append(t)
        for r in range(self_dimensions[0]):
            for c in range(self_dimensions[1]):
                temp_vector[r][c] = self.vector[r][c] * -1

        return Matrix(temp_vector, field=self.field)

    def __rmul__(self, second_matrix):
        if type(second_matrix) == int or type(second_matrix) == float:
            self_dimensions = self.get_dimension()
            temp_vector = []
            for r in range(self_dimensions[0]):
                t = []
                for c in range(self_dimensions[1]):
                    t.append(0)
                temp_vector.append(t)
            for r in range(self_dimensions[0]):
                for c in range(self_dimensions[1]):
                    temp_vector[r][c] = self.vector[r][c] * second_matrix

            return Matrix(temp_vector, field=self.field)
        else:
            raise Exception("Unsupported operation")

    def __mul__(self, second_matrix):
        if type(second_matrix) == int or type(second_matrix) == float:
            self_dimensions = self.get_dimension()
            temp_vector = []
            for r in range(self_dimensions[0]):
                t = []
                for c in range(self_dimensions[1]):
                    t.append(0)
                temp_vector.append(t)
            for r in range(self_dimensions[0]):
                for c in range(self_dimensions[1]):
                    temp_vector[r][c] = self.vector[r][c] * second_matrix

            return Matrix(temp_vector, field=self.field)
        elif type(second_matrix) == Matrix:
            if self.field == second_matrix.field:
                self_dimensions = self.get_dimension()
                second_matrix_dimension = second_matrix.get_dimension()

                if self_dimensions[1] == second_matrix_dimension[0]:
                    temp_vector = []
                    for r in range(self_dimensions[0]):
                        t = []
                        for c in range(second_matrix_dimension[1]):
                            t.append(0)
                        temp_vector.append(t)
                    
                    for r in range(self_dimensions[0]):
                        for c in range(second_matrix_dimension[1]):
                            value = self.field(0)

                            for n in range(self_dimensions[1]):
                                value += (self.vector[r][n]*second_matrix.vector[n][c])

                            temp_vector[r][c] = value
                
                    return Matrix(temp_vector, field=self.field)
                else:
                    raise Exception(f"Cannot multiply matrix of dimension {self_dimensions} with matrix of dimension {second_matrix_dimension}")
                    
            else:
                raise Exception("They are not in the same domain")
        else:
            raise Exception("Unsupported operation")

    def change_field(self, new_field):
        if new_field == None:
            raise Exception("The new field cannot be None")
        
        rows, columns = self.get_dimension()
        return_vector = deepcopy(self.vector)

        for r in range(rows):
            for c in range(columns):
                return_vector[r][c] = new_field(return_vector[r][c].value)
        
        return Matrix(return_vector, field=new_field)
        
    def prepend_row(self, r):
        if type(r) == Matrix:
            r = r.vector[0]

        rows, columns = self.get_dimension()
        if len(r) != columns:
            raise Exception(f"Cannot prepend a row of different column size: {len(r)} != {columns}")

        return_vector = []
        return_vector.append(r)

        for l in self.vector:
            return_vector.append(l)

        return Matrix(return_vector, field=self.field)

    def prepend_column(self, c):
        T = self.transpose()
        T = T.prepend_row(c)
        T = T.transpose()

        return T
        
    def append_row(self, r):
        if type(r) == Matrix:
            r = r.vector[0]
        rows, columns = self.get_dimension()
        if len(r) != columns:
            raise Exception(f"Cannot append a row of different column size: {len(r)} != {columns}")

        return_vector = deepcopy(self.vector)
        return_vector.append(r)

        return Matrix(return_vector, field=self.field)

    def append_column(self, c):
        T = self.transpose()
        T = T.append_row(c)
        T = T.transpose()

        return T

    def l_1_norm(self):
        rows, columns = self.get_dimension()
        if columns == 1:
            sum = 0
            for r in range(rows):
                sum += self.vector[r][0]
            return sum.value
        else:
            raise Exception("You can call this method only on vectors (dim: n x 1)")

    def l_2_norm(self):
        rows, columns = self.get_dimension()
        if columns == 1:
            sum = 0
            for r in range(rows):
                sum += (self.vector[r][0]*self.vector[r][0])
            return sqrt(sum.value)
        else:
            raise Exception("You can call this method only on vectors (dim: n x 1)")


    def l_inf_norm(self):
        rows, columns = self.get_dimension()
        if columns == 1:
            max = 0
            for r in range(rows):
                if abs(self.vector[r][0].value) > max:
                    max = abs(self.vector[r][0].value)
            return max
        else:
            raise Exception("You can call this method only on vectors (dim: n x 1)")

    def __getitem__(self, key):
        row, column = key
        if row is None and column is None:
            raise Exception("You must specify the indexes")
        elif column == None:
            return Matrix(self.vector[row], self.field)
        elif row == None:
            rows = self.get_dimension()[0]
            N = Random_Matrix_Generator.zero_matrix(rows, 1)
            for r in range(rows):
                N[r,0] = self[r,column]
            N.field = self.field
            return N
        else:
            return self.vector[row][column]

    def __mod__(self, modulo):
        new_field = Zmod(modulo)
        return_vector = deepcopy(self.vector)

        return Matrix(return_vector, field=new_field)

    def __setitem__(self, key, value):
        row, column = key
        self.vector[row][column] = self.field(value)

    # This function returns the minor matrix of the element
    # at postion i,j
    def get_minor_matrix(self,i,j):
        minor_matrix = [row[:j] + row[j+1:] for row in (self.vector[:i]+self.vector[i+1:])]
        return Matrix(minor_matrix, field=self.field)


    def real(self):
        return_vector = deepcopy(self.vector)
        integers = ZZ()
        return Matrix(return_vector, field=integers)

    def element_wise_product(self, M):
        if self.get_dimension() == M.get_dimension() and self.field == M.field:
            return_vector = deepcopy(self.vector)
            rows, columns = self.get_dimension()
            for r in range(rows):
                for c in range(columns):
                    return_vector[r][c] *= M[r,c]

            return Matrix(return_vector, self.field)
        else:
            raise Exception("Different rows and columns between the two matrices or different fields")


    # This function find the determinant of the matrix
    def determinant(self):
        #base case for 2x2 matrix
        rows, columns = self.get_dimension()
        if rows != columns:
            raise Exception("you cannot compute the determinant of a non square matrix")
        if rows == 2 and columns == 2:
            return self.vector[0][0]*self.vector[1][1]-self.vector[0][1]*self.vector[1][0]
        determinant = 0
        for c in range(rows):
            determinant += ((-1)**c)*self.vector[0][c]*(self.get_minor_matrix(0,c)).determinant()
        return determinant

    def getMatrixInverse(self):
        rows, columns = self.get_dimension()

        if rows != columns:
            raise Exception("you cannot compute the inverse of a non square matrix")
        
        deter = self.determinant()
        if(deter==0):
            raise Exception('Error! Determinant of the matrix is zero.\nInverse cannot be calculated.')
        
        #special case for 2x2 matrix:
        if rows == 2 and columns == 2:
            return [[self.vector[1][1]/deter, -1*self.vector[0][1]/deter],
                    [-1*self.vector[1][0]/deter, self.vector[0][0]/deter]]
        #find matrix of cofactors
        cofactors = []
        for r in range(rows):
            cofactorRow = []
            for c in range(columns):
                minor = self.get_minor_matrix(r,c)
                cofactorRow.append(((-1)**(r+c)) * minor.determinant())
            cofactors.append(cofactorRow)
            
        # Transposing the cofactor matrix to get disjoint
        cofactor_matrix = Matrix(cofactors, self.field)
        cofactor_matrix = cofactor_matrix.transpose()
        
        # Dividing each element of the adjoint matrix by the determinant
        cofactor_rows, cofactor_cols = cofactor_matrix.get_dimension()
        for r in range(cofactor_rows):
            for c in range(cofactor_cols):
                cofactor_matrix.vector[r][c] = cofactor_matrix.vector[r][c]/deter
                
        return cofactor_matrix

    def __eq__(self, M2):
        if self.get_dimension() != M2.get_dimension():
            return False
        
        if self.field != M2.field:
            return False

        rows, columns = self.get_dimension()
        for r in range(rows):
            for c in range(columns):
                if self[r,c].value != M2[r,c].value:
                    return False
        return True

        


class Random_Matrix_Generator:
    def from_uniform_distribution(B: int, rows: int, columns: int, field):
        matrix = []
        for _ in range(rows):
            t = [field(random.randint(-(B//2), B//2)) for i in range(columns)]
            matrix.append(t)
        return Matrix(matrix, field)

    def from_bernoulli_distribution(B: int, rows: int, columns: int, field):
        matrix = []
        for _ in range(rows):
            t = [field(sample_bernoulli_exp(Fraction(0.5), random)) for i in range(columns)]
            matrix.append(t)
        return Matrix(matrix, field)

    def from_gaussian_distribution(B: int, mu, sigma, rows: int, columns: int, field):
        matrix = []
        for _ in range(rows):
            #t = [field(int(random.gauss(mu,sigma))) for i in range(columns)]
            t = [field(int(sample_dgauss(sigma))) for i in range(columns)]

            matrix.append(t)
        return Matrix(matrix, field=field)        

    def zero_matrix(rows: int, columns: int):
        matrix = []
        for _ in range(rows):
            t = []
            for _ in range(columns):
                t.append(0)
            matrix.append(t)
        
        field=ZZ()
        return Matrix(matrix, field=field)

    def identity_matrix(n: int, field):
        matrix = []
        for i in range(n):
            t = []
            for j in range(n):
                if i == j:
                    t.append(1)
                else:
                    t.append(0)
            matrix.append(t)
        
        return Matrix(matrix, field=field)
    






