#!/usr/bin/sage -python

from sage.all import GF, prod, EllipticCurve, PolynomialRing

class TestECMont:

    def test_xADD_formula(self):
        """Test that the formula used for xADD works correctly"""

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

        XP = P.x()
        XQ = Q.x() 
        XR = R.x()
        ZP = ZQ = ZR = 1

        # xADD formula in projective coordinates: P + Q = (X' : Z')
        # X' = ZR * (XP * XQ - ZP * ZQ)^2 
        # Z' = XR * (XP * ZQ - XQ * ZP)^2
        X = ZR * (XP * XQ - ZP * ZQ)**2
        Z = XR * (XP * ZQ - XQ * ZP)**2

        x = X/Z
        assert x == T.x()

        # Test in projective space
        # Multiply by 'random' scalars => should still work
        XP, ZP = [ (416*i + 7)   * coord for coord in [XP, ZP] ]
        XQ, ZQ = [ (247*i + 159) * coord for coord in [XQ, ZQ] ]
        XR, ZR = [ (75*i + 79)   * coord for coord in [XR, ZR] ]

        X = ZR * (XP * XQ - ZP * ZQ)**2
        Z = XR * (XP * ZQ - XQ * ZP)**2
        x = X/Z
        assert x == T.x()

    def test_xISOG_formula(self):
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
