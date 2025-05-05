#!/usr/bin/sage

from sage.all import EllipticCurve, Primes, randint, order_from_multiple, gcd, is_prime, prod, GF
from verifiers.protocols import MSIDH 

def verify_test_msidh_gen_pub_params():
    print("---: test_msidh_gen_pub_params()")

    t = 4
    print("t:", t)
    p, A, B, f = MSIDH.gen_pub_params(t)
    print("p:", p)
    print("A:", A)
    print("B:", B)
    print("f:", f)
    assert A == 4 * 5
    assert B == 3 * 7
    assert f == 1
    assert p == A * B * f - 1

    t = 9
    print("t:", t)
    p, A, B, f = MSIDH.gen_pub_params(t)
    print("p:", p)
    print("A:", A)
    print("B:", B)
    print("f:", f)
    assert A == 4 * 5 * 11 * 17 * 23
    assert B == 3 * 7 * 13 * 19
    assert f == 2
    assert p == A * B * f - 1

class TestT4:

    @classmethod
    def load_globals(cls):
        """Utility function to use in sage environment after 'load' is called"""
        global i, p, F, E, A, B, f
        i = cls.i
        p = cls.p
        F = cls.F
        E = cls.E
        A = cls.A 
        B = cls.B
        f = cls.f

    @classmethod
    def setup_params(cls):
        """Prepare variables for running tests"""
        # i is used for constructing fp2 elements: a + b * i

        # p = 419, A = 20, B = 21, f = 1
        cls.p, cls.A, cls.B, cls.f = MSIDH.gen_pub_params(4)
        cls.n = cls.p + 1

        assert cls.p == 419
        assert cls.A == 20
        assert cls.B == 21
        assert cls.f == 1

        cls.F = GF(cls.p ** 2, modulus=[1,0,1], names="i")
        # Unpack the tuple using comma after i -> ',', same as [0] syntax
        cls.i, = cls.F._first_ngens(1)

        # Montgomery Starting Curve E: y^2 = x^3 + 6x^2 + x
        cls.E = EllipticCurve(cls.F, [0, 6, 0, 1, 0])

        if not cls.E.is_supersingular():
            raise ValueError("E is not a supersingular curve")

    @classmethod
    def verify_test_msidh_gen_pubkey(cls):
        print("---: test_msidh_gen_pubkey()")
        i, E, A, B, n = cls.i, cls.E, cls.A, cls.B, cls.n

        P = E(295*i + 398, 158*i + 219)
        Q = E(314*i + 149, 240*i + 165)
        PQd = E(29*i + 395, 58*i + 150)
        assert PQd == P - Q

        print(f"xP: {P.x()}")
        print(f"xQ: {Q.x()}")
        print(f"xPQd: {PQd.x()}")

        # Check if valid full torsion basis
        assert n == 420
        assert order_from_multiple(P.weil_pairing(Q, n), n, operation="*") == n

        # Construct Alice basis (PA, QA) = [n//A](P, Q)
        PA = E(128*i + 414, 248*i + 385)
        QA = E(269*i + 32, 396*i + 61)
        PQAd = E(107*i + 300, 387*i + 41)

        print(f"xPA: {PA.x()}")
        print(f"xQA: {QA.x()}")
        print(f"xPQAd: {PQAd.x()}")
        assert PA == P * (n // A)
        assert QA == Q * (n // A)
        assert PA - QA == (P - Q) * (n // A) == PQAd
        assert order_from_multiple(PA.weil_pairing(QA, A), A, operation="*") == A

        # Construct Bob basis (PB, QB) = [n//B](P, Q)
        PB = E(214*i + 177, 156*i + 347)
        QB = E(65*i + 173, 79*i + 203)
        PQBd = E(173*i + 197, 405*i + 98)

        print(f"xPB: {PB.x()}")
        print(f"xQB: {QB.x()}")
        print(f"xPQBd: {PQBd.x()}")
        assert PB == P * (n // B)
        assert QB == Q * (n // B)
        assert PB - QB == (P - Q) * (n // B) == PQBd
        assert order_from_multiple(PB.weil_pairing(QB, B), B, operation="*") == B

        # 1, 8, 13 20: are the quadratic roots of unity modulo 21 (B)
        for j, (a_sec, a_mask) in enumerate([(2, 1), (13, 8), (3, 13), (15, 20)]):
            assert pow(a_mask, 2, B) == 1
            print(f"a_sec: {a_sec}")
            print(f"a_mask: {a_mask}")

            msidh = MSIDH(cls.p, A, B, E, PA, QA, a_sec, a_mask, is_bob=False, mont_model=True)
            AE, APB, AQB = msidh.gen_pubkey(PB, QB)
            a2_24p = (AE.a2() + 2)/4
            print(f"aφ(E)(24p): {a2_24p}")
            print(f"xφ(PB): {APB.x()}")
            print(f"xφ(QB): {AQB.x()}")
            print(f"xφ(PQBd): {(APB-AQB).x()}")


if __name__ == "__main__":

    verify_test_msidh_gen_pub_params()

    # Small MSIDH t = 4
    TestT4.setup_params()
    TestT4.verify_test_msidh_gen_pubkey()

