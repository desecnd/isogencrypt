
from sage.all import factor, order_from_multiple, randrange, CRT_list, EllipticCurve, prod 
def validate_torsion_basis(P, Q, n: int) -> bool:
    """Return True if points (P, Q) create valid (non-degenerate) torsion basis of order n"""
    return order_from_multiple(P.weil_pairing(Q, n), n, operation='*') == n

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

def _mont_coef_2(K):
    """Return 2-isogeny codomain coefficient A in Montgomery Model
    
    Raises:
        ValueError: Given kernel is point (0,0) of order 2 leading to singular curve.
    """

    if K.order() != 2 or K.x() == 0:
        raise ValueError("Kernel is invalid - will land on singular curve")

    return 2 * (1 - 2 * K.x() ** 2)

def _mont_coef_odd(K):
    """Return Montgomery A coefficient of the odd-degree isogenous codomain curve given kernel K"""
    assert K.order() % 2 == 1

    d = (K.order() - 1)//2
    kpts = [ K * j for j in range(1, d + 1) ]

    # A' is the coefficient of the codomain
    A = K.curve().a2()
    sigma = sum([KP.x() for KP in kpts])
    sigma_hat = sum([1/KP.x() for KP in kpts])
    pi = prod(KP.x() for KP in kpts)
    A_ = (6 * (sigma_hat - sigma) + A) * pi**2

    return A_

def mont_coef(K):
    """Return montgomery curve coefficient Ax^2 given kernel of the isogeny of any degree"""

    # If point at inf. return curent curve coeff
    if K == K.curve()(0):
        return K.curve().a2()

    mul = K.order()
    K0 = K

    F = K.curve().base()
    E_ord = (F.characteristic() + 1) ** 2

    # Prime factors with repetition (i.e: 2^2 -> [2, 2])
    prime_factors = [ p for (p, m) in factor(mul) for _ in range(m) ]

    # For each prime degree of the K order construct isogeny to montgomery model curve
    for p in prime_factors:
        EK = K0.curve()

        # Use formulas for 2-isogeny and odd-degree isogeny
        T = K0 * (mul // p)
        a = _mont_coef_2(T) if p == 2 else _mont_coef_odd(T)

        # Construct the codomain curve with calculated coefficient
        ET = EllipticCurve(F, [0, a, 0, 1, 0])
        ET.set_order(E_ord)

        # Push the point to the next curve
        K0 = EK.isogeny(T, codomain=ET, algorithm="traditional")(K0)
        mul //= p

    return a

def mont_isog(K):
    """Construct isogeny with codomain equal to the montgomery x-only arithmetic (internally using Weierstrass model)"""
    E0 = K.curve()
    A = mont_coef(K)
    E1 = EllipticCurve(E0.base(), [0, A, 0, 1, 0])
    return E0.isogeny(K, codomain=E1)