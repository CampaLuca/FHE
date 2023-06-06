import sys
sys.path.insert(0, '../')

from GSW.Messages_in_small_space.GSW import *

#########################################
#                 TEST                  #
#########################################

# security parameter (best value: 24)
k = 24

# gaussian distribution (you can change them and see what are the differences)
mu    = 0.0
sigma = 1

DEBUG = True
# gsw initialization
gsw = GSW(k, mu, sigma, debug=DEBUG)
gsw.key_gen()


print("Running the tests\n")
print("Status: ", end="", flush=True)

counter = 0
for _ in range(100):
    check = True

    m1 = random.randint(0,1)
    m2 = random.randint(0,1)

    C1 = gsw.encrypt(m1)
    C2 = gsw.encrypt(m2)
    C3 = gsw.homomorphic_addition(C1, C2)
    C4 = gsw.homomorphic_multiplication(C1, C2)
    C5 = gsw.NAND(C1, C2)

    check = check and gsw.decrypt(C1) == m1
    check = check and gsw.decrypt(C2) == m2

    check = check and gsw.decrypt(C3) == (m1+m2)
    check = check and gsw.decrypt(C4) == m1*m2
    check = check and gsw.decrypt(C5) == 1 - m1*m2

    if check:
        counter += 1

    print("#", end="", flush=True)


print(f"Correct tests: {counter}/100")


