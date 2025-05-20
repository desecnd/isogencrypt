#!/usr/bin/sage -python

from sage.all import EllipticCurve, GF, prod
from isogencrypt_sage.isogeny import mont_coef, mont_isog, _mont_coef_2, _mont_coef_odd

class TestP431:

    @classmethod
    def load_globals(cls):
        """Load global variables used often in SageMath and other tests to omit calling to cls"""

        # Bind from class to global variables
        global p, F, i, E
        p, F, i, E = cls.p, cls.F, cls.i, cls.E 

    @classmethod
    def setup_class(cls):
        """Prepare global variables for running tests"""
        # p + 1 = 2^4 * 3^3
        p = 431
        F = GF(p ** 2, modulus=[1,0,1], names="i")
        i, = F._first_ngens(1)

        # Montgomery Starting Curve E: y^2 = x^3 + 6x^2 + x
        E = EllipticCurve(F, [0, 6, 0, 1, 0])
        assert E.is_supersingular()

        # Bind to global variables
        cls.p, cls.F, cls.i, cls.E = p, F, i, E

    @classmethod
    def test_xDBL_small(cls):
        print("---: test_xDBL_small")
        cls.load_globals()

        # Point with order 432
        P = E(292 + 15 * i, 281 + 235 * i)
        Q = E(61 + 184 * i, 395 + 90 * i)

        # Q = [2]P
        assert 2 * P == Q
        assert Q.x() == 61 + 184 * i

    @classmethod
    def test_xADD_small(cls):
        print("---: test_xADD_small")
        cls.load_globals()

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

    @classmethod
    def test_criss_cross_small(cls):
        print("---: test_criss_cross_small")
        cls.load_globals()

        x = F(416*i + 175)
        y = F(112*i + 179)
        z = F(235*i + 107)
        w = F(183*i + 197)

        e = x * w + y * z
        f = x * w - y * z
        
        assert e == 367*i + 314
        assert f == 19*i + 425

        print(f"x: {x}")
        print(f"y: {y}")
        print(f"z: {z}")
        print(f"w: {w}")
        print(f"xw+yz: {e}")
        print(f"xw-yz: {f}")

        
    @classmethod
    def test_xLADDER3PT(cls):
        print("---: test_xLADDER3PT")
        cls.load_globals()

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

    @classmethod
    def test_point_normalize_coords(cls):
        print("---: test_point_normalize_coords()")
        cls.load_globals()

        X = 395*i + 201
        print(f"X: {X}")
        Y = 272*i + 286
        print(f"Y: {Y}")

        X_ = X / Y
        print(f"X': {X_}")
        assert X_ == 12*i + 95

    @classmethod
    def test_xDBLe(cls):
        print("---: test_xDBLe")
        cls.load_globals()

        P = E(387*i + 387, 325*i + 125)
        assert P.order() == 432
        print(f"xP: {P.x()}")

        P2 = E(400*i + 311, 412*i + 256)
        print(f"x[2]P: {P2.x()}")
        assert P2 == P * 2

        P4 = E(13*i + 67, 206*i + 377)
        print(f"x[4]P: {P4.x()}")
        assert P4 == P * 4

        P8 = E(213*i + 105, 373*i + 392)
        print(f"x[8]P: {P8.x()}")
        assert P8 == P * 8

        P2_12345 = E(304*i + 223, 134*i + 11)
        print(f"x[2^1235]P: {P2_12345.x()}")
        assert P2_12345 == P * (2 ** 12345)


    @classmethod
    def test_ISOG2e(cls):
        print("---: test_ISOG2e")
        cls.load_globals()

        K0 = E(33*i + 429, 205*i + 374)
        K = K0
        assert K.order() == 2**4
        print(f"xK: {K.x()}")

        # Push-through point
        P0 = E(158*i + 183, 43*i + 20)
        P = P0
        print(f"xP: {P.x()}")
        # Has full torsion and lays above different 2-order point than K
        assert P.order() == 432
        assert P * (432//2) != K * (16//2)

        Ei = E
        for j in range(4):
            # T = [2^(e - 1)]K
            T = (K * 2 ** (4 - j - 1)) if j < 3 else K
            A = mont_coef(T)
            E_ = EllipticCurve(F, [0, A, 0, 1, 0])
            phi = Ei.isogeny(T, codomain=E_)
            P = phi(P)
            Ei = E_
            if j < 3:
                K = phi(K)

        # Additional check that the formula for mont_coef works
        # Only x-coordinate is taken into account in mont_isog fomula, to it can differ
        P_ = mont_isog(K0)(P0)
        assert P == P_ or P == -P_ 
        assert A == mont_coef(K0)

        print(f"xφ(P): {P.x()}")
        print(f"aφ(E): {A}")

class TestP139:

    @classmethod
    def load_globals(cls):
        # Bind from class to global variables
        global p, F, i, E
        p, F, i, E = cls.p, cls.F, cls.i, cls.E 

    @classmethod
    def setup_class(cls):
        """Prepare global variables for running tests"""
        # p + 1 = 2^2 * 5 * 7 
        p = 139

        F = GF(p ** 2, modulus=[1,0,1], names="i")
        i, = F._first_ngens(1)

        # Montgomery Starting Curve E: y^2 = x^3 + 6x^2 + x
        E = EllipticCurve(F, [0, 6, 0, 1, 0])
        assert E.is_supersingular()

        # Bind to global variables
        cls.p, cls.F, cls.i, cls.E = p, F, i, E


    @classmethod
    def test_KPS(cls):
        print("---: test_KPS()")
        cls.load_globals()

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


    @classmethod
    def test_xISOG_and_aISOG(cls): 
        print("---: test_xISOG_and_aISOG()")
        cls.load_globals()

        # Define isogeny kernel of order 5
        K = E(77*i + 38, 87*i + 133)
        print(f"xK: {K.x()}")
        assert K.order() == 5

        P = E(32*i + 42, 97*i + 88)
        # Full-order point
        assert P.order() == 140
        print(f"xP: {P.x()}")

        # Montgomery Curve with A_ coefficient
        A_ = _mont_coef_odd(K)
        print(f"aφ(E): {A_}")

        E_ = EllipticCurve(F, [0, A_, 0, 1, 0])
        assert E_.j_invariant() == 100

        phi = E.isogeny(K, codomain=E_, algorithm='traditional')
        phi_P = phi(P)
        print(f"xφ(P): {phi_P.x()}")

        # Order of the point did not change
        # assert phi_P.order() == 140

    @classmethod
    def test_xLADDER(cls): 
        print("---: test_xLADDER()")
        cls.load_globals()

        P = E(7*i + 97, 9*i + 129) 
        print(f"xP: {P.x()}")

        P1 = E(98*i + 43, 23*i + 67)
        print(f"x[2^80]P: {P1.x()}")
        assert P1 == P * 2**80

        P2 = E(56*i + 96, 14*i + 39)
        print(f"x[2^80-1]P: {P2.x()}")
        assert P2 == P * (2**80 - 1)

        P3 = E(94*i + 31, 80*i + 23)
        m = 0xf5697b000f01c17d4c5e
        print(f"x[{hex(m)}]P: {P3.x()}")
        assert P3 == P * m

    @classmethod
    def test_xLADDER_int(cls): 
        print("---: test_xLADDER_int()")
        cls.load_globals()

        P = E(108*i + 136, 68*i + 134)
        print(f"xP: {P.x()}")

        muls = [ (P * j).x() for j in range(1, 5)] 
        assert muls == [108*i + 136, 113*i + 131, 42*i + 83, 47*i + 107]

        for m, mulx in enumerate(muls):
            print(f"x[{m + 1}]P = {mulx}")

    @classmethod
    def test_ISOG_chain_odd(cls):
        print("---: test_ISOG_chain_odd()")
        cls.load_globals()

        K = E(108*i + 136, 68*i + 134)
        assert K.order() == 35
        print(f"xK: {K.x()}")

        K5 = E(96*i + 71, 15*i + 87) 
        assert K5.order() == 5
        assert K5 == K * 7
        print(f"xK5: {K5.x()}")

        # K5.isogeny()
        A5 = _mont_coef_odd(K5)
        E5 = EllipticCurve(F, [0, A5, 0, 1, 0])
        print(f"aE5: {A5}")

        phi5 = E.isogeny(K5, codomain=E5)
        K7 = phi5(K)
        print(f"xφ(K): {K7.x()}")

        assert K7.order() == 7

        A7 = _mont_coef_odd(K7)
        E7 = EllipticCurve(F, [0, A7, 0, 1, 0])
        phi7 = E5.isogeny(K7, codomain=E7)

        assert phi7(K7).is_zero()
        print(f"aE7: {A7}")

    @classmethod
    def test_xISOG2_and_aISOG2(cls):
        print("---: test_xISOG2_and_aISOG2()")
        cls.load_globals()

        K2 = E(100*i + 136, 0)
        assert K2.order() == 2
        print(f"xK2: {K2.x()}")

        P = E(70*i + 36, 129*i + 11)
        assert P.order() == 140
        print(f"xP: {P.x()}")

        # Formula for obtaining the A' for 2-isogeny
        A2 = 2 * (1 - 2 * K2.x() ** 2)
        print(f"aE2: {A2}")
        A2_24p = (A2 + 2) / 4
        print(f"aE2(24p): {A2_24p}")

        E2 = EllipticCurve(F, [0, A2, 0, 1, 0])
        phi2 = E.isogeny(K2, codomain=E2)
        
        P_ = phi2(P)
        print(f"xφ(P): {P_.x()}")

    @classmethod
    def test_ISOG2_bad_point_error(cls):
        print("---: test_ISOG2_bad_point_error()")
        cls.load_globals()

        K2 = E(0, 0)
        assert K2.order() == 2
        print(f"xK: {K2.x()}")

        # Formula for obtaining the A' for 2-isogeny
        # This will give value 2 which is incorrect (singular curve)
        A2 = 2 * (1 - 2 * K2.x() ** 2)
        print(f"aφ(K): {A2}")

        try:
            EllipticCurve(F, [0, A2, 0, 1, 0])
            assert False and "Constructed curve should be singular"
        except ArithmeticError:
            pass

    @classmethod
    def test_ISOG_chain(cls):
        print("---: test_ISOG_chain()")
        cls.load_globals()

        K0 = E(34*i + 99, 25*i + 95) 
        print(f"xK: {K0.x()}")
        assert K0.order() == 140
        # K140 does not lie over (0, 0) point
        assert (70 * K0).x() != 0

        P0 = E(8*i + 137, 51*i + 35)
        print("xP:", P0.x())
        assert P0.order() == 140
        assert (70 * P0).x() != 0
        # P0 lays over different order 2 point than K0
        assert (70 * P0) != (70 * K0)

        # 1. Calculate phi1: E -> E1 with ker = T1 of degree 2 
        T1 = K0 * 70
        assert T1.order() == 2

        # a) Calculate codomain
        A1 = _mont_coef_2(T1)

        # b) Construct the isogeny
        E1 = EllipticCurve(F, [0, A1, 0, 1, 0])
        phi1 = E.isogeny(T1, codomain=E1)

        # c) Push the kernel K through the isogeny
        K1 = phi1(K0)
        assert K1.order() == 70
        P = phi1(P0)

        # 2. Calculate phi2: E1 -> E2 with ker = T2 of degree 2 
        T2 = K1 * 35
        assert T2.order() == 2

        # a) Calculate codomain
        A2 = _mont_coef_2(T2)

        # b) Construct the isogeny
        E2 = EllipticCurve(F, [0, A2, 0, 1, 0])
        phi2 = E1.isogeny(T2, codomain=E2)

        # c) Push the kernel K through the isogeny
        K2 = phi2(K1)
        assert K2.order() == 35
        P = phi2(P)

        # 3. Calculate phi3: E2 -> E3 with ker = T3 of degree 5
        T3 = K2 * 7
        assert T3.order() == 5

        # a) Calculate codomain
        A3 = _mont_coef_odd(T3)

        # b) Construct the isogeny
        E3 = EllipticCurve(F, [0, A3, 0, 1, 0])
        phi3 = E2.isogeny(T3, codomain=E3)

        # c) Push the kernel K through the isogeny
        K3 = phi3(K2)
        assert K3.order() == 7
        P = phi3(P)

        # 4. Calculate phi4: E3 -> E4 with ker = T4 of degree 7
        T4 = K3
        assert T4.order() == 7

        # a) Calculate codomain
        A4 = _mont_coef_odd(T4)

        # b) Construct the isogeny
        E4 = EllipticCurve(F, [0, A4, 0, 1, 0])
        phi4 = E3.isogeny(T4, codomain=E4)

        # c) Push the kernel K through the isogeny
        K4 = phi4(K3)
        assert K4 == E4(0)

        P = phi4(P)
        A = A4

        # Additional check that the formula for mont_coef works
        # Only x-coordinate is taken into account in mont_isog fomula, to it can differ
        P_ = mont_isog(K0)(P0)
        assert P == P_ or P == -P_ 
        assert A == mont_coef(K0)

        print(f"aφ(E): {A}")
        print(f"xφ(P): {P.x()}")

    @classmethod
    def test_j_invariant(cls):
        print("---: test_j_invariant()")
        cls.load_globals()

        A = F(92*i + 25)
        j = F(9*i + 29)
        assert j == EllipticCurve(F, [0, A, 0, 1, 0]).j_invariant()
        assert j == 256 * (A**2 - 3)**3 / (A**2 - 4)
        print(f"a(E): {A}")
        print(f"j(E): {j}")

        A = F(125*i + 99)
        j = F(79*i + 30)
        assert j == EllipticCurve(F, [0, A, 0, 1, 0]).j_invariant()
        assert j == 256 * (A**2 - 3)**3 / (A**2 - 4)
        print(f"a(E): {A}")
        print(f"j(E): {j}")

        A = F(43*i + 61)
        j = F(78*i + 97)
        assert j == EllipticCurve(F, [0, A, 0, 1, 0]).j_invariant()
        assert j == 256 * (A**2 - 3)**3 / (A**2 - 4)
        print(f"a(E): {A}")
        print(f"j(E): {j}")



def main():
    # SIDH-like prime: 
    TestP431.setup_class()
    TestP431.test_xDBL_small()
    TestP431.test_xADD_small()
    TestP431.test_criss_cross_small()
    TestP431.test_xLADDER3PT()
    TestP431.test_point_normalize_coords()
    TestP431.test_xDBLe()
    TestP431.test_ISOG2e()

    # Odd-degree prime:
    TestP139.setup_class()
    TestP139.test_KPS()
    TestP139.test_xISOG_and_aISOG()
    TestP139.test_xLADDER_int()
    TestP139.test_xLADDER()
    TestP139.test_ISOG_chain_odd()
    TestP139.test_xISOG2_and_aISOG2()
    TestP139.test_ISOG2_bad_point_error()
    TestP139.test_ISOG_chain()
    TestP139.test_j_invariant()

if __name__ == '__main__':
    main()
