#!/usr/bin/sage

from math import prod
from sage.all import *

# Method of torsion basis sampling based on one by grhkm
def sample_torsion_basis_smooth(E, r):
    assert E.is_supersingular()

    # p + 1 for supersingular curve
    p = E.base().characteristic()
    assert (p + 1)**2 == E.order()

    # Point order must divide the E order
    assert (p + 1) % r == 0

    # for now we assume that r is a product of primes
    # TODO: think how to make this example work for exp > 1
    assert all(e == 1 for _, e in factor(r))

    # we don't care about the case where  
    coeff = (p + 1) / r
    # gcd(coeff, r) # this does not matter at all

    def _sample_torsion_point():
        full_torsion = False
        while not full_torsion:

            # Get point that order is divided by at least 1 prime divisor li
            while True:
                X = E.random_point() * coeff
                if X != E(0):
                    break

            # Check for each prime divisor li if it divides the point order
            for li, _ in factor(r):
                if X * (r // li) == E(0):
                    break
            else:
                # Point has full r-torsion for composite r
                full_torsion = True

        return X

    # Start with P = Q
    P = Q = _sample_torsion_point()

    # r is composite, therefore we can by accident sample element which 
    # multiplicative order only divides r and is not r 
    # in such case (P, Q) do not create a basis and isogeny degree
    # will be different than r
    #
    # We are using "order_from_multiple" instead of "e.multiplicative_order()" due
    # to the high inefficiency of the latter in SageMath 10.4
    while (e := P.weil_pairing(Q, r)) == 1 or order_from_multiple(e, r, operation='*') != r:
        # If degenerate case (same subgroup), repeat sampling
        Q = _sample_torsion_point()

    return P, Q

def sample_quadratic_root_of_unity(modulus: int):
    """
    Sample an element x from Z/bZ where x ** 2 = 1 mod b
    """

    # TODO: I guess this works only for the list of primes with exp = 1 max? 
    # If prime^2 divides n then something might not work?

    # For each factor, find the square roots
    roots = []
    moduli = []

    # Assume factorizatoin is easy
    for prime, exp in factor(modulus):
        assert exp == 1 and "[!] What if exp > 1?"

        # 1 or p - 1 -> Toss a coin
        root = prime - 1 if randrange(2) == 0 else 1

        roots.append(root)
        moduli.append(prime)
            
    # Use CRT to find the root
    res = CRT_list(roots, moduli)
    assert (res ** 2) % modulus == 1

    return res


# t-bit security
t = 64
P = Primes()
ll = [ P.unrank(i) for i in range(1, 2 * t + 1, 2) ]
qq = [ P.unrank(i) for i in range(2, 2 * t + 1, 2) ]

assert len(ll) == t
assert len(qq) == t

A = prod(ll)
B = prod(qq)
AB = A * B

assert gcd(A, B) == 1

for f in range(1000):
    p = f * AB - 1
    if is_prime(p):
        print(f"{f = }")
        break
else:
    raise ValueError("Cannot find cofactor for p")

# Find Supersingular Elliptic Curve over Fp^2
F.<i> = GF(p**2)
E0 = EllipticCurve(j=supersingular_j(F))
E0.set_order((p+1)**2)
assert E0.is_supersingular()

print(f"[*] Calculating Torsion Basis...(PA, QA)")
PA, QA = sample_torsion_basis_smooth(E0, A)
assert PA.order() == QA.order() == A 
assert order_from_multiple(PA.weil_pairing(QA, A), A, operation='*') == A

print(f"[*] Calculating Torsion Basis...(PB, QB)")
PB, QB = sample_torsion_basis_smooth(E0, B)
assert PB.order() == QB.order() == B
# This is highly inefficient
# assert PB.weil_pairing(QB, B).multiplicative_order() == B
assert order_from_multiple(PB.weil_pairing(QB, B), B, operation='*') == B

# By defuault weil_pairing result is a element_pari_ffelt type
# which order returns additive order

# Calculate Alice Kernel
a = randrange(2, A)
A_ker = PA + a * QA
assert A_ker.order() == A

# Calculate Bob Kernel
b = randrange(2, B)
B_ker = PB + b * QB
assert B_ker.order() == B

# Calculate isogeny for Alice 
A_phi = E0.isogeny(A_ker, algorithm='factored')

# Calculate isogeny for Bob
B_phi = E0.isogeny(B_ker, algorithm='factored')

print(f"{A_phi.degree() / A = }")
print(f"{B_phi.degree() / B = }")

# Calculate public key of Alice
# Sample alpha for modulus B  
alpha = sample_quadratic_root_of_unity(B)
EA = A_phi.codomain()

# Send basis points (PB, QB) through the isogeny
A_phiPB = alpha * A_phi(PB)
A_phiQB = alpha * A_phi(QB)
A_pub = (EA, A_phiPB, A_phiQB)

# Calculate public key of Bob
beta = sample_quadratic_root_of_unity(A)
EB = B_phi.codomain()

B_phiPA = beta * B_phi(PA)
B_phiQA = beta * B_phi(QA)
B_pub = (EB, B_phiPA, B_phiQA)

# Calculate Alice shared key
# 1. Perform the pairing check for basis

print("[*] Calculating Alice shared secret")

# Isogeny of coprime degree does not change the order of the points
assert B_phiPA.order() == B_phiQA.order() == A

ea = B_phiPA.weil_pairing(B_phiQA, A) 
if ea == PA.weil_pairing(QA, A) ** B:
    print("[+] Bob pk valid")
else:
    print("[!] Bob pk invalid")

A_ker2 = B_phiPA + a * B_phiQA
A_phi2 = EB.isogeny(A_ker2, algorithm='factored')

EBA = A_phi2.codomain()
A_ss = EBA.j_invariant()

print("[*] Calculating Bob shared secret")
# Isogeny of coprime degree does not change the order of the points
assert B_phiPA.order() == B_phiQA.order() == A

eb = A_phiPB.weil_pairing(A_phiQB, B) 
if eb == PB.weil_pairing(QB, B) ** A:
    print("[+] Alice pk valid")
else:
    print("[!] Alice pk invalid")

B_ker2 = A_phiPB + b * A_phiQB
B_phi2 = EA.isogeny(B_ker2, algorithm='factored')

EAB = B_phi2.codomain()
B_ss = EAB.j_invariant()

assert B_ss == A_ss
print(f"{A_ss = }")
print(f"{B_ss = }")
