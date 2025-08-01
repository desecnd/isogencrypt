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
    def test_point_normalize_coords(cls):
        print("---: test_point_normalize_coords()")
        cls.load_globals()

        X = 395*i + 201
        print(f"X: {X}")
        Y = 272*i + 286
        print(f"Y: {Y}")

        x = X / Y
        print(f"x: {x}")
        assert x == 12*i + 95

    @classmethod
    def test_xDBL_small(cls):
        print("---: test_xDBL_small()")
        cls.load_globals()

        # Point with order 432
        P = E(292 + 15 * i, 281 + 235 * i)
        Q = E(61 + 184 * i, 395 + 90 * i)
        print(f"xP: {P.x()}")

        assert 2 * P == Q
        print(f"x[2]P: {Q.x()}")
        
    @classmethod
    def test_xDBLe(cls):
        print("---: test_xDBLe()")
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
        print(f"x[2^12345]P: {P2_12345.x()}")
        assert P2_12345 == P * (2 ** 12345)

    @classmethod
    def test_xADD_small(cls):
        print("---: test_xADD_small()")
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
    def test_xLADDER3PT(cls):
        print("---: test_xLADDER3PT()")
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
            print(f"x[{m + 1}]P: {mulx}")

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
    TestP431.test_point_normalize_coords()
    TestP431.test_xDBL_small()
    TestP431.test_xDBLe()
    TestP431.test_xADD_small()
    TestP431.test_xLADDER3PT()

    # Odd-degree prime:
    TestP139.setup_class()
    TestP139.test_xLADDER_int()
    TestP139.test_xLADDER()
    TestP139.test_j_invariant()

if __name__ == '__main__':
    main()
