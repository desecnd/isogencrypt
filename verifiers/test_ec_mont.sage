#!/usr/bin/sage

from sage.all import GF, prod, is_prime, EllipticCurve, PolynomialRing

def find_ell_prime_number(factors: list[int]) -> int:
    """Find prime number that contains specific factors under of p + 1 and makes a valid supersingular ec graph"""

    ITERS = 100
    for i in range(1, ITERS):
        p = prod(factors) * 2 * i - 1
        if not is_prime(p):
            continue
        try:
            # Modulus must be x^2 + 1, otherwise the arithmetic may be wrong
            U = GF(p**2, modulus=[1,0,1], names="u")
            u, = U._first_ngens(1)
            if (EllipticCurve(U, [1, 0]).is_supersingular()):
                return p
        except:
            continue
    else:
        raise ValueError(f"Cannot find valid prime number in {ITERS} iterations.")

def test_xADD_formula():
    """Test that the formula used for xADD works correctly"""

    print("TEST ---: xADD")
    p = 431
    F = GF(p**2, modulus=[1,0,1], names="i")
    i, = F._first_ngens(1)
    # Montgomery starting curve E: y^2 = x^3 + 6x^2 + x
    E = EllipticCurve(F, [0, 6, 0, 1, 0])

    P = E(107*i + 11, 301*i + 428)
    Q = E(71*i + 243, 61*i + 18)

    R = P - Q 
    assert R == E(286*i + 119, 301*i + 124)

    T = P + Q
    assert T == E(62*i + 260, 42*i + 180)

    print(f"x(P + Q) = {T.x()}")

    print("1) (X : Z) = (x : 1)")
    XP = P.x()
    XQ = Q.x() 
    XR = R.x()
    ZP = ZQ = ZR = 1

    # xADD formula in projective coordinates: P + Q = (X' : Z')
    # X' = ZR * (XP * XQ - ZP * ZQ)^2 
    # Z' = XR * (XP * ZQ - XQ * ZP)^2
    X = ZR * (XP * XQ - ZP * ZQ)**2
    Z = XR * (XP * ZQ - XQ * ZP)**2

    print(f"{X = }")
    print(f"{Z = }")
    x = X/Z
    print(f"> X/Z = {x = }")
    assert x == T.x()

    print("2) (X : Z) = (ax : a)")
    # Test in projective space
    # Multiply by 'random' scalars => should still work
    XP, ZP = [ (416*i + 7)   * coord for coord in [XP, ZP] ]
    XQ, ZQ = [ (247*i + 159) * coord for coord in [XQ, ZQ] ]
    XR, ZR = [ (75*i + 79)   * coord for coord in [XR, ZR] ]

    X = ZR * (XP * XQ - ZP * ZQ)**2
    Z = XR * (XP * ZQ - XQ * ZP)**2
    print(f"{X = }")
    print(f"{Z = }")
    x = X/Z
    print(f"> X/Z = {x = }")
    assert x == T.x()

def test_xISOG_formula():
    print("TEST ---: xISOG")
    # Test if we obtain the exact same formula for f(x) polynomial 
    # given isogeny of degree 5

    # Define the parameters
    # p + 1 = 4 * 5 * 7
    p = 139
    F = GF(p ** 2, modulus=[1,0,1], names="i")
    i, = F._first_ngens(1)
    A = 6
    # E: y^2 = x^3 + Ax^2 + x 
    E = EllipticCurve(F, [0, A, 0, 1, 0])
    assert E.is_supersingular()

    # K is the kernel of the 5-degree isogeny
    K = E(77*i + 38, 87*i + 133)
    assert K.order() == 5

    # Formulas from: https://eprint.iacr.org/2017/504.pdf
    # Equation (5):
    sigma = sum([ (K * j).x() for j in range(1, 3) ])
    sigma_inv = sum([ 1/(K * j).x() for j in range(1, 3) ])
    pi = prod([ (K * j).x() for j in range(1, 3) ])

    # A' is the coefficient of the codomain
    A_ = (6 * sigma_inv - 6 * sigma + A) * pi**2
    # B' does not matter when only x-arithmetic is used
    # B_ = B * pi^2

    # Codomain elliptic curve E' calculated from formulas above
    E_ = EllipticCurve(F, [0, A_, 0, 1, 0])

    # Define polynomial ring with cooeficients in Fp^2
    # Equation (6):
    R = PolynomialRing(F, names="x")
    x, = R._first_ngens(1)
    x_map = x * prod([ ((x * (K * j).x() - 1)/(x - (K * j).x()))**2 for j in range(1, 3) ])
    print(f"x_map: {x_map}")

    # Use Velu formulas to obtain isogeny of degree 5
    phi = E.isogeny(K, algorithm='traditional')

    # We arrrive at the same curve under j-inv but different coefficients
    assert phi.codomain().j_invariant() == E_.j_invariant()
    assert phi.codomain().a_invariants() != E_.a_invariants()

    # Isogeny as a rational map has the same denominator as <K> is the same
    assert phi.x_rational_map().denominator() == x_map.denominator()
    assert phi.x_rational_map().numerator() != x_map.numerator()

    # We must create isogeny
    psi = phi.codomain().isomorphism_to(E_)

    # FInal isogeny to E' is equal to Montgomery-form formulas
    iso = psi * phi
    assert iso.codomain() == E_
    assert iso.x_rational_map() == x_map

    # We can tell Sage to calculate isogeny explicitly to E'
    # Velu formula, but joined with isomorphism to E'
    assert iso == E.isogeny(K, codomain=E_)

if __name__ == '__main__':
    test_xADD_formula()
    test_xISOG_formula()