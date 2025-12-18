#!/usr/bin/sage -python

from sage.all import EllipticCurve, GF
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
    def test_criss_cross_small(cls):
        print("---: test_criss_cross_small()")
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
    def test_ISOG2e(cls):
        print("---: test_ISOG2e()")
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

        print(f"deg: {K.order()}")
        print(f"n: {K.order()//2}")

        # Montgomery Curve with A_ coefficient
        A_ = _mont_coef_odd(K)
        print(f"aφ(K): {A_}")

        E_ = EllipticCurve(F, [0, A_, 0, 1, 0])
        assert E_.j_invariant() == 100

        phi = E.isogeny(K, codomain=E_, algorithm='traditional')
        phi_P = phi(P)
        print(f"xφ(P): {phi_P.x()}")

        # Order of the point did not change
        assert phi_P.order() == 140

    @classmethod
    def test_ISOG_chain_odd(cls):
        print("---: test_ISOG_chain_odd()")
        cls.load_globals()

        K = E(108*i + 136, 68*i + 134)
        print(f"xK: {K.x()}")
        assert K.order() == 35

        K5 = E(96*i + 71, 15*i + 87) 
        # print(f"xK5: {K5.xy()}")
        assert K5.order() == 5
        assert K5 == K * 7

        # Construct next montgomery curve from calculated codomain coefficient
        A5 = _mont_coef_odd(K5)
        E5 = EllipticCurve(F, [0, A5, 0, 1, 0])
        assert A5 == 76*i + 85

        phi5 = E.isogeny(K5, codomain=E5)
        K7 = phi5(K)
        # print(f"xφ(K): {K7.xy()}")
        assert K7 == E5((30*i + 60, 125*i + 136))
        assert K7.order() == 7

        A7 = _mont_coef_odd(K7)
        E7 = EllipticCurve(F, [0, A7, 0, 1, 0])
        phi7 = E5.isogeny(K7, codomain=E7)

        assert phi7(K7).is_zero()
        print(f"aφ(K): {A7}")

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

        print(f"aφ(K): {A}")
        print(f"xφ(P): {P.x()}")

    @classmethod
    def test_ISOG_chain_trivial(cls):
        print("---: test_ISOG_chain_trivial()")
        cls.load_globals()

        K = E(0)
        print(f"XK: {K[0]}")
        print(f"ZK: {K[0]}")

        phi = E.isogeny(K)
        assert phi.degree() == 1
        assert phi.codomain() == E

        P0 = E(8*i + 137, 51*i + 35)
        print("xP:", P0.x())
        Q0 = phi(P0)

        a_coeff = phi.codomain().a_invariants()[1]
        print(f"aφ(K): {a_coeff}")
        print(f"xφ(P): {Q0.x()}")



def main():
    # SIDH-like prime: 
    TestP431.setup_class()
    TestP431.test_criss_cross_small()
    TestP431.test_ISOG2e()

    # Odd-degree prime:
    TestP139.setup_class()
    TestP139.test_KPS()
    TestP139.test_xISOG_and_aISOG()
    TestP139.test_ISOG_chain_odd()
    TestP139.test_xISOG2_and_aISOG2()
    TestP139.test_ISOG_chain()
    TestP139.test_ISOG_chain_trivial()

if __name__ == '__main__':
    main()
