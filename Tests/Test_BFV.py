import sys
sys.path.insert(0, '../')


from BFV.BFV import *

#########################################
#                 TEST                  #
#########################################

### bfv params
p = 32  # p should be a power of 2
N = 1024 # the degree of the polynomials will be N-1

# gaussian distribution
mu    = 0
sigma = 0.5 * 3.2

# bfv initialization
bfv = BFV(N, p, mu, sigma)
bfv.key_gen()



print("Running the tests\n")
print("Status: ", end="", flush=True)
counter = 0
for _ in range(100):
    check = True
    ###### Creating two messages
    m1, m2 = randint(0,2**(p-1)-1), randint(0,2**(p-1)-1)

    M_1 = binary_encode_message(m1, N)
    M_2 = binary_encode_message(m2, N)
    M1 = Polynomial_ring(p, Polynomial(N, M_1), bfv.ntt_context) 
    M2 = Polynomial_ring(p, Polynomial(N, M_2), bfv.ntt_context) 

    #print(f'M1 = {m1}, M2 = {m2}')

    # Checking encryption and decryption
    m1_after_decryption = decode_message(bfv.decrypt(bfv.encrypt(M1)).poly.coefficients, p)
    check = check and m1 == m1_after_decryption


    # Checking homomorphic addition
    C1 = bfv.encrypt(M1)
    C2 = bfv.encrypt(M2)
    C3 = bfv.homomorphic_addition(C1, C2)
    sum = m1 + m2
    #print(f'M1 + M2 = {sum}')
    sum_after_homomorphic_addition = bfv.decrypt(C3)
    check = check and sum == decode_message(sum_after_homomorphic_addition.poly.coefficients, p)

    # Checking naive homomorphic multiplication (with SK^2)
    multiplication = m1*m2
    #print(f'M1 * M2 = {multiplication}')
    C4 = bfv.naive_homomorphic_multiplication(C1, C2)
    multiplication_after_naive_homomorphic_mult = bfv.decrypt(C4)
    check = check and multiplication == decode_message(multiplication_after_naive_homomorphic_mult.poly.coefficients, p)

    # Checking homomorphic multiplication with RELINEARIZATION, without BASE DECOMPOSITION
    # it does not work because the value C_3 * e increases the noise making the decryption infeasible
    """C5 = bfv.homomorphic_multiplication(C1, C2, base_decomposition=False)
    multiplication_after_relinearization_homomorphic_mult = bfv.decrypt(C5)
    assert multiplication == decode_message(multiplication_after_relinearization_homomorphic_mult.poly.coefficients, p)"""

    # Checking homomorphic multiplication with RELINEARIZATION AND BASE DECOMPOSITION
    C6 = bfv.homomorphic_multiplication(C1, C2, base_decomposition=True)
    multiplication_after_relinearization_base_decomp_homomorphic_mult = bfv.decrypt(C6)
    check = check and multiplication == decode_message(multiplication_after_relinearization_base_decomp_homomorphic_mult.poly.coefficients, p)

    if check:
        counter += 1

    print("#", end="", flush=True)


print()
print(f"Correct tests: {counter}/100")
