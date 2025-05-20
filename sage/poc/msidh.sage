#!/usr/bin/sage
#
# PoC MSIDH implementation in Weierstrass Curve Model in SageMath

from sage.all import EllipticCurve, factor, Primes, randrange, order_from_multiple, CRT_list, gcd, is_prime, GF, prod

def sample_torsion_basis_smooth(E, r: int, montgomery_basis: bool = False):
    """Fast method for finding smooth torsion basis on supersingular elliptic curve"""

    assert E.is_supersingular()

    # Grab characteristic from the base field E is defined over
    p = E.base().characteristic()
    assert (p + 1)**2 == E.order()

    # Point order must divide the E order
    assert (p + 1) % r == 0

    # Prime factors with repetition (i.e: 2^2 -> [2, 2])
    prime_factors = [ p for (p, m) in factor(r) for _ in range(m) ]

    # We don't care about the case where gcd(coeff, r) != 1 -> this does not matter at all
    coeff = (p + 1) / r

    def _sample_torsion_point():
        while True: 

            # Get point that order is divided by at least 1 prime divisor li
            while True:
                X = E.random_point() * coeff
                if X != E(0):
                    break

            # Check for each prime divisor li if it divides the point order
            for li in prime_factors:
                if X * (r // li) == E(0):
                    break
            else:
                # Point has full r-torsion for composite r -> break
                break

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
    while True:
        Q = _sample_torsion_point()

        e = P.weil_pairing(Q, r)
        if e == 1:
            continue

        e_ord = order_from_multiple(e, r, operation='*')
        if e_ord != r:
            continue

        # Montgomery 2-isogeny requires that Q lays above point (0, 0)
        if montgomery_basis and r % 2 == 0:
            # Q lays above the point (0, 0), we can break
            if (Q * (r // 2)).x() == 0:
                break

            # P lays above the point (0, ), we can swap the points
            if (P * (r // 2)).x() == 0:
                P, Q = Q, P
                break

            # We need to search further
            continue
    
        # No "wrong" conditions -> we can break
        break

    return P, Q

def sample_quadratic_root_of_unity(modulus: int):
    """
    Sample an element x from Z/bZ where x ** 2 = 1 mod b
    """

    # For each factor, find the square roots
    roots = []
    moduli = []

    # Assume factorizatoin is easy
    for prime, exp in factor(modulus):
        # This can work for n = 4, so n does not have to be prime
        n = prime ** exp

        # 1 or n - 1 -> Toss a coin
        root = n - 1 if randrange(2) == 0 else 1

        roots.append(root)
        moduli.append(n)
            
    # Use CRT to find the root
    res = CRT_list(roots, moduli)
    assert (res ** 2) % modulus == 1

    return res


if __name__ == '__main__':
    # Number of primes to multiply  ~ t/2 bit classical and ~t/4 quantum security
    t = 30
    print(f"[%] Running MSIDH with t={t}")

    P = Primes(proof=False)
    ll = [ P.unrank(2 * i) for i in range((t+1)//2) ]
    qq = [ P.unrank(2 * i + 1) for i in range(t//2) ]
    # Alice uses 4 as the first number instead of 2 (probably due to quadratic root -1 = 1 for 2, but not for 4)
    ll[0] = 4
    assert len(ll) + len(qq) == t

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

    # Same as: F.<i> = GF(p**2)
    F = GF(p**2, names=('i',), modulus=[1, 0, 1])
    (i,) = F._first_ngens(1)

    E0 = EllipticCurve(F, [0, 6, 0, 1, 0])
    E0.set_order((p+1)**2)
    assert E0.is_supersingular()

    print("[*] Calculating Torsion Basis...(PA, QA)")
    PA, QA = sample_torsion_basis_smooth(E0, A)
    assert PA.order() == QA.order() == A 
    assert order_from_multiple(PA.weil_pairing(QA, A), A, operation='*') == A

    print("[*] Calculating Torsion Basis...(PB, QB)")
    PB, QB = sample_torsion_basis_smooth(E0, B)
    assert PB.order() == QB.order() == B
    # This is highly inefficient
    # assert PB.weil_pairing(QB, B).multiplicative_order() == B
    assert order_from_multiple(PB.weil_pairing(QB, B), B, operation='*') == B

    # By default weil_pairing result is a element_pari_ffelt type
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
