from sage.all import factor, order_from_multiple, randrange, CRT_list, EllipticCurve, prod

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

    # we don't care about the case where gcd(coeff, r) != 1 -> this does not matter at all
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
    """Return 2-isogeny codomain coefficient A in Montgomery Model"""
    assert K.order() == 2 and K.x() != 0
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
    """Combine odd and even montgomery coeff functions into general formula"""
    n = K.order()
    K_ = K 

    # Extract even isogeny by combining 2-isogenies
    while n % 2 == 0:
        EK = K_.curve()
        T = K_ * (n//2)
        A = _mont_coef_2(T)
        ET = EllipticCurve(EK.base(), [0, A, 0, 1, 0])
        K_ = EK.isogeny(T, codomain=ET)(K_)
        n //= 2

    # Based on the calculated kernel, solve odd_isogeny
    return _mont_coef_odd(K_)