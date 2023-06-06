import sys
sys.path.insert(0, '../')

from GSW.Messages_in_bigger_space.GSW import *

#########################################
#                 TEST                  #
#########################################

### bfv params
k = 24  # the number of bits

# gaussian distribution
mu    = 0.0
sigma = 1

DEBUG = True
# bfv initialization
gsw = GSW(k, mu, sigma, debug=DEBUG)
gsw.key_gen()


print("Running the tests\n")
print("Status: ", end="", flush=True)

counter = 0
for _ in range(100):
    check = True

    m1 = random.randint(0,k)
    m2 = random.randint(0,k)
    sum = m1 + m2
    mult = m1*m2
    C1 = gsw.encrypt(m1)
    C2 = gsw.encrypt(m2)
    check = check and  gsw.decrypt(C1) == m1
    check = check and  gsw.decrypt(C2) == m2
    check = check and  gsw.decrypt(gsw.homomorphic_addition(C1, C2)) == sum
    check = check and  gsw.decrypt(gsw.homomorphic_multiplication(C1, C2)) == mult

    if check:
        counter += 1
    
    print("#", end="", flush=True)


print(f"Correct tests: {counter}/100")
