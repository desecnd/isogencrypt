#!/usr/bin/sage 

# Sage-friendly definition in global scope 
# otherwise i is bound to the class and defining 
# numbers "x + y*i" with cls.i or self.i gets tedious
p = None
F = None 
i = None 
E = None

def _montgomery_coef(K):
    """Return Montgomery A coefficient of the odd-degree isogenous codomain curve given kernel K"""
    assert K.order() % 2 == 1

    d = (K.order() - 1)//2
    kpts = [ K * j for j in range(1, d + 1) ]

    # A' is the coefficient of the codomain
    A = K.curve().a2()
    sigma = sum([KP.x() for KP in kpts])
    sigma_hat = sum([1/KP.x() for KP in kpts])
    pi = prod(KP.x() for KP in kpts)
    A_ = (6 * sigma_hat - 6 * sigma + A) * pi^2

    return A_

class TestcaseP431:

    def setup_params():
        """Prepare global variables for running tests"""
        global p, F, i, E

        # p + 1 = 2^4 * 3^3
        p = 431
        F.<i> = GF(p^2, modulus=[1,0,1])
        # Montgomery Starting Curve E: y^2 = x^3 + 6x^2 + x
        E = EllipticCurve(F, [0, 6, 0, 1, 0])

        if not E.is_supersingular():
            raise ValueError("E is not a supersingular curve")

    def verify_test_xDBL_small():
        print("VERIFY ---: test_xDBL_small")
        # Point with order 432
        P = E((292 + 15 * i, 281 + 235 * i))
        Q = E(61 + 184 * i, 395 + 90 * i)

        # Q = [2]P
        assert 2 * P == Q
        assert Q.x() == 61 + 184 * i

    # Function for testing xADD:
    def verify_test_xADD_small():
        print("VERIFY ---: test_xADD_small")

        P = E(271*i + 259, 422*i + 97)
        print(f"xP: {P.x()}")

        Q = E(335*i + 262, 69*i + 198)
        print(f"xQ: {Q.x()}")

        PQdiff = E(411*i + 143, 245*i + 213)
        assert PQdiff == P - Q

        print(f"xP-Q: {PQdiff.x()}")
        PQsum = P + Q

        assert PQsum == E(106*i + 416 , 111*i + 405)
        print(f"xP+Q: {PQsum.x()}")

        assert PQsum.x() == 416 + 106 * i

    def verify_test_criss_cross_small():
        print("VERIFY ---: test_criss_cross_small")
        x = F(416*i + 175)
        y = F(112*i + 179)
        z = F(235*i + 107)
        w = F(183*i + 197)

        e = x * w + y * z
        f = x * w - y * z
        # assert e == 206*i + 256
        # assert f == 180*i + 52 
        print(x * w)
        print(y * z)

        print(f"x: {x}")
        print(f"y: {y}")
        print(f"z: {z}")
        print(f"w: {w}")
        print(f"xw+yz: {e}")
        print(f"xw-yz: {f}")

        
    def verify_test_xLADDER3PT_small():
        print("VERIFY ---: test_xLADDER3PT_small")
        P = E(271*i + 259, 422*i + 97)
        print(f"xP: {P.x()}")
        Q = E(335*i + 262, 69*i + 198)
        print(f"xQ: {Q.x()}")
        PQd = P - Q
        print(f"xP-Q: {PQd.x()}")
        n = 87
        print(f"n: {n}")

        R = P + n * Q
        print(f"x(P+nQ): {R.x()}")
        assert R == E(45*i + 360, 249*i + 429)

        T = P - n * Q
        print(f"x(P-nQ): {T.x()}")

    def verify_test_point_normalize_coords():
        print("VERIFY ---: test_point_normalize_coords()")

        X = 395*i + 201
        print(f"X: {X}")
        Y = 272*i + 286
        print(f"Y: {Y}")

        X_ = X / Y
        print(f"X': {X_}")
        assert X_ == 12*i + 95

class TestcaseP139:

    def setup_params():
        """Prepare global variables for running tests"""
        global p, F, i, E

        # p + 1 = 2^2 * 5 * 7 
        p = 139
        F.<i> = GF(p^2, modulus=[1,0,1])
        # E = EllipticCurve(F, [1, 0])
        E = EllipticCurve(F, [0, 6, 0, 1, 0])

        if not E.is_supersingular():
            raise ValueError("E is not a supersingular curve")

    def verify_test_KPS():
        print("VERIFY ---: test_KPS()")
        # Define isogeny kernel of order 7
        K = E(101*i + 20, 102*i + 21)
        print(f"xK: {K.x()}")
        assert K.order() == 7

        d = (7 - 1)//2

        # Multiplies of the kernel point
        x_coords = [ (K * j).x() for j in range(1, d + 1) ]
        print(f"xK1: {x_coords[0]}")
        print(f"xK2: {x_coords[1]}")
        print(f"xK3: {x_coords[2]}")
        assert x_coords == [101*i + 20, 82*i + 16, 106*i + 124]


    def verify_test_xISOG_and_aISOG(): 
        print("VERIFY ---: test_xISOG_and_aISOG()")
        # Define isogeny kernel of order 5
        K = E(77*i + 38, 87*i + 133)
        print(f"xK: {K.x()}")
        assert K.order() == 5

        P = E(32*i + 42, 97*i + 88)
        # Full-order point
        assert P.order() == 140
        print(f"xP: {P.x()}")

        # Montgomery Curve with A_ coefficient
        A_ = _montgomery_coef(K)
        print(f"aφ(E): {A_}")

        E_ = EllipticCurve(F, [0, A_, 0, 1, 0])
        assert E_.j_invariant() == 100

        phi = E.isogeny(K, codomain=E_, algorithm='traditional')
        phi_P = phi(P)
        print(f"xφ(P): {phi_P.x()}")

        # Order of the point did not change
        # assert phi_P.order() == 140

    def verify_test_xLADDER(): 
        print("VERIFY ---: test_xLADDER()")

        P = E(108*i + 136, 68*i + 134)
        print(f"xP: {P.x()}")

        muls = [ (P * j).x() for j in range(1, 5)] 
        assert muls == [108*i + 136, 113*i + 131, 42*i + 83, 47*i + 107]

        for m, mulx in enumerate(muls):
            print(f"x[{m + 1}]P = {mulx}")



    def verify_test_ISOG_chain():
        print("VERIFY ---: test_ISOG_chain()")
        K = E(108*i + 136, 68*i + 134)
        assert K.order() == 35
        print(f"xK: {K.x()}")

        K5 = E(96*i + 71, 15*i + 87) 
        assert K5.order() == 5
        assert K5 == K * 7
        print(f"xK5: {K5.x()}")

        # K5.isogeny()
        A5 = _montgomery_coef(K5)
        E5 = EllipticCurve(F, [0, A5, 0, 1, 0])
        print(f"aE5: {A5}")

        phi5 = E.isogeny(K5, codomain=E5)
        K7 = phi5(K)
        print(f"xφ(K): {K7.x()}")

        assert K7.order() == 7

        A7 = _montgomery_coef(K7)
        E7 = EllipticCurve(F, [0, A7, 0, 1, 0])
        phi7 = E5.isogeny(K7, codomain=E7)

        assert phi7(K7).is_zero()
        print(f"aE7: {A7}")

    def verify_test_xISOG2_and_aISOG2():
        print("VERIFY ---: test_xISOG2_and_aISOG2()")

        K2 = E(100*i + 136, 0)
        assert K2.order() == 2
        print(f"xK2: {K2.x()}")

        P = E(70*i + 36, 129*i + 11)
        assert P.order() == 140
        print(f"xP: {P.x()}")

        # Formula for obtaining the A' for 2-isogeny
        A2 = 2 * (1 - 2 * K2.x() ** 2)
        print(f"aE2: {A2}")
        E2 = EllipticCurve(F, [0, A2, 0, 1, 0])
        phi2 = E.isogeny(K2, codomain=E2)
        
        P_ = phi2(P)
        print(f"xφ(P): {P_.x()}")

if __name__ == '__main__':

    # SIDH-like prime: 
    # TestcaseP431.setup_params()
    # TestcaseP431.verify_test_xDBL_small()
    # TestcaseP431.verify_test_xADD_small()
    # TestcaseP431.verify_test_criss_cross_small()
    # TestcaseP431.verify_test_xLADDER3PT_small()
    # TestcaseP431.verify_test_point_normalize_coords()

    # Odd-degree prime:
    TestcaseP139.setup_params()
    TestcaseP139.verify_test_KPS()
    TestcaseP139.verify_test_xISOG_and_aISOG()
    TestcaseP139.verify_test_xLADDER()
    TestcaseP139.verify_test_ISOG_chain()
    TestcaseP139.verify_test_xISOG2_and_aISOG2()

    # TestcaseP139.verify_test_xISOG4()
    # TestcaseP139.verify_test_aISOG4()





