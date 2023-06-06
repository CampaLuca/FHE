import sympy
from MathObj.Number.NTT import modinv


"""
INTEGERS
"""
class ZZ:
    def __init__(self):
        self.domain = 'ZZ'

    def __call__(self, value):
        if type(value) == ZZ_instance or type(value) == Zmod_instance:
            return ZZ_instance(self, value.value)
        elif type(value) == QQ_instance:
            return ZZ_instance(self, round(value.value, 0))
        else:
            return ZZ_instance(self, value)


    def __eq__(self, field):
        if type(field) == ZZ:
            return True
        else:
            return False

class ZZ_instance:
    def __init__(self, field, value):
        self.field = field
        self.value = value

    def __add__(self, value):

        if type(value) == int:
            value = self.field(value)
            return self.field(self.value + value.value)
        elif type(value) == float:
            qq_field = QQ()
            value = qq_field(value)
            return qq_field(self.value + value.value)
        elif type(value) == ZZ_instance:
            return self.field(self.value + value.value)
        elif type(value) == QQ_instance:
            qq_field = QQ()
            return qq_field(self.value + value.value)
        else:
            raise Exception("Operation not supported")

    def __radd__(self, value):
        
        if type(value) == int:
            value = self.field(value)
            return self.field(self.value + value.value)
        elif type(value) == float:
            qq_field = QQ()
            value = qq_field(value)
            return qq_field(self.value + value.value)
        elif type(value) == ZZ_instance:
            return self.field(self.value + value.value)
        elif type(value) == QQ_instance:
            qq_field = QQ()
            return qq_field(self.value + value.value)
        else:
            raise Exception("Operation not supported")

    def __sub__(self, value):
        
        if type(value) == int:
            value = self.field(value)
            return self.field(self.value - value.value)
        elif type(value) == float:
            qq_field = QQ()
            value = qq_field(value)
            return qq_field(self.value - value.value)
        elif type(value) == ZZ_instance:
            return self.field(self.value - value.value)
        elif type(value) == QQ_instance:
            qq_field = QQ()
            return qq_field(self.value - value.value)
        else:
            raise Exception("Operation not supported")
    

    def __rsub__(self, value):
        
        if type(value) == int:
            value = self.field(value)
            return self.field(self.value - value.value)
        elif type(value) == float:
            qq_field = QQ()
            value = qq_field(value)
            return qq_field(self.value - value.value)
        elif type(value) == ZZ_instance:
            return self.field(self.value - value.value)
        elif type(value) == QQ_instance:
            qq_field = QQ()
            return qq_field(self.value - value.value)
        else:
            raise Exception("Operation not supported")

    
    def __mul__(self, value):
        
        if type(value) == int:
            value = self.field(value)
            return self.field(self.value * value.value)
        elif type(value) == float:
            qq_field = QQ()
            value = qq_field(value)
            return qq_field(self.value * value.value)
        elif type(value) == ZZ_instance:
            return self.field(self.value * value.value)
        elif type(value) == QQ_instance:
            qq_field = QQ()
            return qq_field(self.value * value.value)
        else:
            raise Exception("Operation not supported")
    
    def __rmul__(self, value):
        
        if type(value) == int:
            value = self.field(value)
            return self.field(self.value * value.value)
        elif type(value) == float:
            qq_field = QQ()
            value = qq_field(value)
            return qq_field(self.value * value.value)
        elif type(value) == ZZ_instance:
            return self.field(self.value * value.value)
        elif type(value) == QQ_instance:
            qq_field = QQ()
            return qq_field(self.value * value.value)
        else:
            raise Exception("Operation not supported")

    def __truediv__(self, value):  
        qq_field = QQ()

        if type(value) == int:
            value = self.field(value)
            return qq_field(self.value / value.value)
        elif type(value) == float:
            value = qq_field(value)
            return qq_field(self.value / value.value)
        elif type(value) == ZZ_instance:
            return qq_field(self.value / value.value)
        elif type(value) == QQ_instance:
            return qq_field(self.value / value.value)
        else:
            raise Exception("Operation not supported")

    def __floordiv__(self, value):
        if type(value) == int:
            value = self.field(value)
            return self.field(self.value // value.value)
        else:
            raise Exception("Operation not supported")

    def __mod__(self, modulo):
        new_field = Zmod(modulo)
        return new_field(self.value)
    

""" INTEGERS MODULO n """
class Zmod:
    def __init__(self, modulo):
        self.modulo = modulo

    def __call__(self, value):
        if type(value) == ZZ_instance or type(value) == Zmod_instance:
            return Zmod_instance(self, value.value % self.modulo)
        elif type(value) == QQ_instance:
            return Zmod_instance(self, round(value.value, 0) % self.modulo)
        else:
            return Zmod_instance(self, value)

    def __eq__(self, field):
        if type(field) == Zmod:
            var_local = vars(self)
            var_field = vars(field)

            for el1 in var_field:
                if var_local[el1] != var_field[el1]:
                    return False
            
            return True
        else:
            return False

class Zmod_instance:
    def __init__(self, mod_field : Zmod, value):
        self.value = value % mod_field.modulo
        self.mod_field = mod_field

    def __add__(self, value):
        if type(value) == int or type(value) == float:
            value = self.mod_field(value)

        if self.mod_field.modulo == value.mod_field.modulo:
            return Zmod_instance(self.mod_field, (self.value + value.value) % self.mod_field.modulo)
        else:
            raise Exception("They are not in the same domain")

    def __radd__(self, value):
        if type(value) == int or type(value) == float:
            value = self.mod_field(value)
            
        if self.mod_field.modulo == value.mod_field.modulo:
            return Zmod_instance(self.mod_field, (self.value + value.value) % self.mod_field.modulo)
        else:
            raise Exception("They are not in the same domain")

    def __sub__(self, value):
        if type(value) == int or type(value) == float:
            value = self.mod_field(value)

        if self.mod_field.modulo == value.mod_field.modulo:
            return Zmod_instance(self.mod_field, (self.value - value.value) % self.mod_field.modulo)
        else:
            raise Exception("They are not in the same domain")

    def __rsub__(self, value):
        if type(value) == int or type(value) == float:
            value = self.mod_field(value)

        if self.mod_field.modulo == value.mod_field.modulo:
            return Zmod_instance(self.mod_field, (self.value - value.value) % self.mod_field.modulo)
        else:
            raise Exception("They are not in the same domain")
    
    def __mul__(self, value):
        if type(value) == int or type(value) == float:
            value = self.mod_field(value)

        if self.mod_field.modulo == value.mod_field.modulo:
            return Zmod_instance(self.mod_field, (self.value * value.value) % self.mod_field.modulo)
        else:
            raise Exception("They are not in the same domain")

    def __rmul__(self, value):
        if type(value) == int or type(value) == float:
            value = self.mod_field(value)

        if self.mod_field.modulo == value.mod_field.modulo:
            return Zmod_instance(self.mod_field, (self.value * value.value) % self.mod_field.modulo)
        else:
            raise Exception("They are not in the same domain")

    def __truediv__(self, value):
        if type(value) == int or type(value) == float:
            value = self.mod_field(value)

        if self.mod_field.modulo == value.mod_field.modulo:
            try:
                inverse_value = modinv(value.value, self.mod_field.modulo)
                return Zmod_instance(self.mod_field, (self.value * inverse_value) % self.mod_field.modulo)
            except:
                raise(f"Modulo inverse of {value.value} modulo {value.mod_field.modulo} does not exist")
        else:
            raise Exception("They are not in the same domain")

    def __mod__(self, modulo):
        new_field = Zmod(modulo)
        return new_field(self.value)

    
""" RATIONAL NUMBERS """
class QQ:
    def __init__(self):
        self.domain = 'QQ'

    def __call__(self, value):
        if type(value) == ZZ_instance or type(value) == Zmod_instance or type(value) == QQ_instance:
            return QQ_instance(self, value.value)
        else:
            return QQ_instance(self, value)

    def __eq__(self, field):
        if type(field) == QQ:
            return True
        else:
            return False

class QQ_instance:
    def __init__(self, field, value):
        self.field = field
        self.value = value

    def __add__(self, value):  
        if type(value) == int or type(value) == float:
            value = self.field(value)
            return self.field(self.value + value.value)
        elif type(value) == QQ_instance or type(value) == ZZ_instance:
            return self.field(self.value + value.value)
        else:
            raise Exception("Operation not supported")

    def __radd__(self, value):
        if type(value) == int or type(value) == float:
            value = self.field(value)
            return self.field(self.value + value.value)
        elif type(value) == QQ_instance or type(value) == ZZ_instance:
            return self.field(self.value + value.value)
        else:
            raise Exception("Operation not supported")

    def __sub__(self, value):
        if type(value) == int or type(value) == float:
            value = self.field(value)
            return self.field(self.value - value.value)
        elif type(value) == QQ_instance or type(value) == ZZ_instance:
            return self.field(self.value - value.value)
        else:
            raise Exception("Operation not supported")
    

    def __rsub__(self, value):
        if type(value) == int or type(value) == float:
            value = self.field(value)
            return self.field(self.value - value.value)
        elif type(value) == QQ_instance or type(value) == ZZ_instance:
            return self.field(self.value - value.value)
        else:
            raise Exception("Operation not supported")

    
    def __mul__(self, value):
        if type(value) == int or type(value) == float:
            value = self.field(value)
            return self.field(self.value * value.value)
        elif type(value) == QQ_instance or type(value) == ZZ_instance:
            return self.field(self.value * value.value)
        else:
            raise Exception("Operation not supported")
    
    def __rmul__(self, value):
        if type(value) == int or type(value) == float:
            value = self.field(value)
            return self.field(self.value * value.value)
        elif type(value) == QQ_instance or type(value) == ZZ_instance:
            return self.field(self.value * value.value)
        else:
            raise Exception("Operation not supported")

    def __truediv__(self, value):  
        if type(value) == int or type(value) == float:
            value = self.field(value)
            return self.field(self.value / value.value)
        elif type(value) == QQ_instance or type(value) == ZZ_instance:
            return self.field(self.value / value.value)
        else:
            raise Exception("Operation not supported")


    