#!/usr/bin/sage

from sage.all import EllipticCurve, order_from_multiple, GF
from isogencrypt_sage.isogeny import validate_torsion_basis 
from isogencrypt_sage.msidh import MSIDH 

def test_msidh_gen_pub_params():
    print(f"---: test_msidh_gen_pub_params()")

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

    t = 20
    print("t:", t)
    p, A, B, f = MSIDH.gen_pub_params(t)
    print("p:", p)
    print("A:", A)
    print("B:", B)
    print("f:", f)

    t = 100
    print("t:", t)
    p, A, B, f = MSIDH.gen_pub_params(t)
    print("p:", p)
    print("A:", A)
    print("B:", B)
    print("f:", f)

    t = 170
    print("t:", t)
    p, A, B, f = MSIDH.gen_pub_params(t)
    print("p:", p)
    print("A:", A)
    print("B:", B)
    print("f:", f)

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
    def setup_class(cls):
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
    def test_msidh_internals(cls):
        print("---: test_msidh_internals()")
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

        # 1, 8, 13 20: are the quadratic roots of unity modulo 21 (B) used by Alice for masking
        # 1, 9, 11, 19: are quadratic roots of unity modulo 20 (A) used by Bob for masking
        for j, (a_sec, a_mask, b_sec, b_mask) in enumerate([(2, 1, 5, 9), (13, 8, 7, 1), (3, 13, 17, 19), (15, 20, 1, 11)]):
            assert pow(a_mask, 2, B) == 1
            assert pow(b_mask, 2, A) == 1

            print(f"a_sec: {a_sec}")
            print(f"a_mask: {a_mask}")
            print(f"b_sec: {b_sec}")
            print(f"b_mask: {b_mask}")

            A_msidh = MSIDH(cls.p, A, B, E, PA, QA, a_sec, a_mask, is_bob=False, mont_model=True)
            AE, APB, AQB = A_msidh.gen_pubkey(PB, QB)
            a2_24p = (AE.a2() + 2)/4
            print(f"aφ(E)(24p): {a2_24p}")
            print(f"xφ(PB): {APB.x()}")
            print(f"xφ(QB): {AQB.x()}")
            print(f"xφ(PQBd): {(APB-AQB).x()}")

            B_msidh = MSIDH(cls.p, A, B, E, PB, QB, b_sec, b_mask, is_bob=True, mont_model=True)
            j_inv = B_msidh.key_exchange(AE, APB, AQB)
            assert B_msidh.tau_ker.order() == B

            EK = B_msidh.EK
            assert j_inv == EK.j_invariant()
            a2_24p = (EK.a2() + 2)/4
            print(f"aτ(φ(E))(24p): {a2_24p}")
            print(f"jτ(φ(E)): {j_inv}")


    @classmethod
    def test_msidh_secret_zero(cls):
        print("---: test_msidh_secret_zero()")
        i, E, n = cls.i, cls.E, cls.n

        P = E(295*i + 398, 158*i + 219)
        Q = E(314*i + 149, 240*i + 165)
        PQd = E(29*i + 395, 58*i + 150)
        assert PQd == P - Q

        print(f"xP: {P.x()}")
        print(f"xQ: {Q.x()}")
        print(f"xPQd: {PQd.x()}")

        # Check if valid full torsion basis
        assert n == 420
        assert validate_torsion_basis(P, Q, n)

        s = 0
        K = P + s * Q
        assert K == P
        print(f"xK: {K.x()}")
        # Print twice to match with the .c test version
        print(f"xK: {K.x()}")

    @classmethod
    def test_msidh_non_deterministic(cls):
        print("---: test_msidh_non_deterministic()")
        i, E, = cls.i, cls.E

        P = E(209*i + 332, 248*i + 396)
        Q = E(345*i + 223, 389*i + 40)
        PQd = E(98*i + 199, 71*i + 181)
        assert PQd == P - Q

        print(f"xP: {P.x()}")
        print(f"xQ: {Q.x()}")
        print(f"xPQd: {PQd.x()}")

    @classmethod
    def test_msidh_monte_carlo(cls):
        print("---: test_msidh_monte_carlo()")
        i, E, = cls.i, cls.E

        P = E(209*i + 332, 248*i + 396)
        Q = E(345*i + 223, 389*i + 40)
        PQd = E(98*i + 199, 71*i + 181)
        assert PQd == P - Q

        print(f"xP: {P.x()}")
        print(f"xQ: {Q.x()}")
        print(f"xPQd: {PQd.x()}")


class TestT30:

    @classmethod
    def load_globals(cls):
        """Utility function to use in sage environment after 'load' is called"""
        global p, A, B, f, n, i, F, E, P, Q
        p, A, B, f, n = cls.p, cls.A, cls.B, cls.f, cls.n 
        i, F, E, P, Q = cls.i, cls.F, cls.E, cls.P, cls.Q 

    @classmethod
    def setup_class(cls):
        # MSIDH params for t = 30
        p, A, B, f = MSIDH.gen_pub_params(30)
        n = p + 1
        assert p == 63220109280835215576290412583087324986549373979
        assert n == A * B * f 
        assert A == 134031250783943894759620
        # assert p == 63156889171554380360714122170504237661562824606019 
        
        a = 6
        F = GF(p ** 2, modulus=[1,0,1], names="i")
        # Unpack the tuple using comma after i -> ',', same as [0] syntax
        i, = F._first_ngens(1)

        # Montgomery Starting Curve E: y^2 = x^3 + 6x^2 + x
        E = EllipticCurve(F, [0, a, 0, 1, 0])
        E.set_order(n**2)

        if not E.is_supersingular():
            raise ValueError("E is not a supersingular curve")

        P = E(32381872305678490404833289490608450363172933700*i + 38566570518230924614310417068523536638018699310, 51454160553705302795561802082188918906573282935*i + 46999345345209075217588309637506871396736048585)
        Q = E(29454235622145096109316297773070819902970029047*i + 17242937661247998401353436361850378272505949076, 4027506100385014783073063257216469725014255007*i + 40503840888595183524703631351014003574434214416)
        P.set_order(n)
        Q.set_order(n)
        # assert validate_torsion_basis(P, Q, n)

        # Q lays over the (0, 0) point of order 2
        assert Q * (n//2) == E(0, 0)

        # Map to class variables
        cls.p, cls.A, cls.B, cls.f, cls.n = p, A, B, f, n
        cls.i, cls.F, cls.E, cls.P, cls.Q = i, F, E, P, Q

    @classmethod
    def test_msidh_internals_large(cls):
        print("---: test_msidh_internals_large()")
        cls.load_globals()

        print(f"xP: {P.x()}")
        print(f"xQ: {Q.x()}")
        print(f"xPQd: {(P - Q).x()}")

        # Construct Alice basis (PA, QA) = [n//A](P, Q)
        PA = P * (n // A)
        QA = Q * (n // A)

        print(f"xPA: {PA.x()}")
        print(f"xQA: {QA.x()}")
        print(f"xPQAd: {(PA - QA).x()}")
        assert validate_torsion_basis(PA, QA, A)

        # Construct Bob basis (PB, QB) = [n//B](P, Q)
        PB = P * (n // B) 
        QB = Q * (n // B)

        print(f"xPB: {PB.x()}")
        print(f"xQB: {QB.x()}")
        print(f"xPQBd: {(PB - QB).x()}")
        assert validate_torsion_basis(PB, QB, B)


        # Secret and Mask for Alice
        a_inputs = [ 
            (11349330453264090324638, 131339333569636892083385),
            (102959793904566459067664, 234614427428290727103469),
        ]

        # Secrets and Masks for Bob
        b_inputs = [
            (72593563569111258510719, 17334875921447289106681),
            (271374038831227351065368, 93531163257802324154959), 
        ]

        for ((A_sec, A_mask), (B_sec, B_mask)) in zip(a_inputs, b_inputs):
            print(f"A_sec: {A_sec}")
            print(f"A_mask: {A_mask}")
            print(f"B_sec: {B_sec}")
            print(f"B_mask: {B_mask}")

            A_msidh = MSIDH(p, A, B, E, PA, QA, A_sec, A_mask, is_bob=False, mont_model=True)
            AE, APB, AQB = A_msidh.gen_pubkey(PB, QB)

            a2_24p = (AE.a2() + 2)/4
            print(f"aφ(E)(24p): {a2_24p}")
            print(f"xφ(PB): {APB.x()}")
            print(f"xφ(QB): {AQB.x()}")
            print(f"xφ(PQBd): {(APB-AQB).x()}")

            B_msidh = MSIDH(p, A, B, E, PB, QB, B_sec, B_mask, is_bob=True, mont_model=True)
            j_inv = B_msidh.key_exchange(AE, APB, AQB)

            EK = B_msidh.EK
            assert j_inv == EK.j_invariant()

            a2_24p = (EK.a2() + 2)/4
            print(f"aτ(φ(E))(24p): {a2_24p}")
            print(f"jτ(φ(E)): {j_inv}")


def main():
    test_msidh_gen_pub_params()

    # Small MSIDH t = 4
    TestT4.setup_class()
    TestT4.test_msidh_internals()
    TestT4.test_msidh_secret_zero()
    TestT4.test_msidh_non_deterministic()
    TestT4.test_msidh_monte_carlo()

    # # Larger MSIDH t = 30
    TestT30.setup_class()
    TestT30.test_msidh_internals_large()

if __name__ == "__main__":
    main()



