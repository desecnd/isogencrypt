#/usr/bin/sage
from sage.all import EllipticCurve, factor, order_from_multiple, Primes, prod, gcd, is_prime, GF,  randint

# Select P and Q such that "none" of them lays above the (0, 0) point
def sample_torsion_basis_smooth(E, r: int, point_above_zero: str = ''):
    """Fast method for finding smooth torsion basis on supersingular elliptic curve
    
    Args:
        point_above_zero: One of ['P', 'Q', 'none'] - determines which point should lay above the (0, 0) point which is important in x-only Montgomery Arithmetic.
        If value other than any of the 3 values is supplied the checks are ignored
    """

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
                # We want to both P and Q be from different torsion subgroup of order 2 than (0, 0)
                if X != E(0):
                    break

            # Check for each prime divisor li if it divides the point order
            for li in prime_factors:
                if X * (r // li) == E(0):
                    break
            else:
                # P cannot lay above the (0, 0) if we supply "none" as argument
                if r % 2 == 0 and point_above_zero == 'none' and (X * (r//2)).x() == 0: 
                    continue

                # Point has full r-torsion for composite r -> break
                break

        return X


    # Start with P = Q
    P = Q = _sample_torsion_point()

    # True if P is in the same order 2 subgrup as E(0, 0)
    P_above_zero = (P * (r // 2)).x() == 0

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
        if r % 2 == 0:
            Q_above_zero = (Q * (r // 2)).x() == 0

            # Select another Q 
            if point_above_zero == 'none' and not Q_above_zero:
                # We should mark this earlier
                assert not P_above_zero
                break

            if P_above_zero or Q_above_zero:
                if point_above_zero == 'P' and not P_above_zero:
                    P, Q = Q, P
                if point_above_zero == 'Q' and not Q_above_zero:
                    P, Q = Q, P

                break

            # We need to search further
            continue
    
        # No "wrong" conditions -> we can break
        break

    return P, Q


def ter_sidh_isogeny(E, PA, QA, primes: list[int], secret: list[int], push_points):
    assert len(secret) == len(primes)
    assert PA.curve() == E
    assert QA.curve() == E

    # P has index 0, 
    # Q has index 1
    cPA, cQA = 1, 1
    nPA, nQA = 1, 1
    for digit, prime in zip(secret, primes):
        if digit == 0:
            nPA *= prime
            cQA *= prime
        elif digit == 1:
            nQA *= prime
            cPA *= prime
        elif digit == 2:
            # Decrease the order of the isogeny
            cPA *= prime
            cQA *= prime
        else:
            raise ValueError("Wrong Digit Value")

    KA1 = PA * cPA
    assert KA1.order() == nPA

    KA2 = QA * cQA
    assert KA2.order() == nQA

    A_phi1 = E.isogeny(KA1, algorithm="factored")
    KA2 = A_phi1(KA2)
    A_phi2 = A_phi1.codomain().isogeny(KA2, algorithm="factored")

    A_phi = A_phi2 * A_phi1

    new_points = [ A_phi(point) for point in push_points ]

    return A_phi.codomain(), new_points


def main():
    t = 30
    print(f"[%] Running terSIDH with t={t}")

    # We will treat "t" differently
    P = Primes(proof=False)
    ll = [ P.unrank(2 * i) for i in range(t) ]
    qq = [ P.unrank(2 * i + 1) for i in range(t) ]
    # Alice uses 4 instead of 2 in order to obtain prime = 3 (mod 4) (p = 4k - 1)
    ll[0] = 4
    assert len(ll) + len(qq) == 2 * t

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


    # Supersingular Elliptic Curve over Fp^2
    # Same as: F.<i> = GF(p**2)
    F = GF(p**2, names=('i',), modulus=[1,0,1])
    (i,) = F._first_ngens(1)

    E0 = EllipticCurve(F, [0, 6, 0, 1, 0])
    E0.set_order((p+1)**2)
    assert E0.is_supersingular()

    print("[*] Calculating Torsion Basis...(PA, QA)")
    PA, QA = sample_torsion_basis_smooth(E0, A, 'none')
    assert PA.order() == QA.order() == A 
    assert order_from_multiple(PA.weil_pairing(QA, A), A, operation='*') == A
    # Points cannot lay above (0, 0) point
    # Otherwise we couldn't use Montgomery Arithmetic, as 
    # 2-isogenies require that kernel does not lay above (0, 0).
    # We are calculating 2-isogenies for each of the P and Q points as kernels
    # This does not matter in Weierstrass Model.
    assert (PA * (A // 2)).x() != 0 
    assert (QA * (A // 2)).x() != 0 

    print("[*] Calculating Torsion Basis...(PB, QB)")
    PB, QB = sample_torsion_basis_smooth(E0, B)
    assert PB.order() == QB.order() == B
    assert order_from_multiple(PB.weil_pairing(QB, B), B, operation='*') == B
    # assert PB * (B // 2) != E0(0, 0) 
    # assert QB * (B // 2) != E0(0, 0)

    # Ternary number
    A_secret = [ randint(0, 2) for _ in range(t) ]
    EA, PQA = ter_sidh_isogeny(E0, PA, QA, ll, A_secret, [PB, QB])
    assert PQA[0].order() == PQA[1].order() == B

    # Apply masking 
    while mA := randint(1, B):
        if gcd(mA, B) == 1:
            break

    APB = mA * PQA[0]
    AQB = pow(mA, -1, B) * PQA[1]

    # B public key calculation
    B_secret = [ randint(0, 2) for _ in range(t) ]
    EB, PQB = ter_sidh_isogeny(E0, PB, QB, qq, B_secret, [PA, QA])
    assert PQB[0].order() == PQB[1].order() == A

    # Apply masking 
    while mB := randint(1, A):
        if gcd(mB, A) == 1:
            break

    BPA = mB * PQB[0]
    BQA = pow(mB, -1, A) * PQB[1]

    # Shared secret calculation
    # Alice:
    EFA, _ = ter_sidh_isogeny(EB, BPA, BQA, ll, A_secret, [])
    # Bob:
    EFB, _ = ter_sidh_isogeny(EA, APB, AQB, qq, B_secret, [])

    assert EFA.j_invariant() == EFB.j_invariant()
    print(EFA.j_invariant())
    print(EFB.j_invariant())


if __name__ == '__main__':
    main()
