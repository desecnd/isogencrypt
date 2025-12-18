#!/usr/bin/sage

from sage.all import EllipticCurve, GF
from isogencrypt_sage.isogeny import validate_torsion_basis 
from isogencrypt_sage.tersidh import TerSIDH

def tersidh_base_valid(P, Q, n: int) -> bool:
    """Return true if (P, Q) is valid torsion basis E[n] for TerSIDH protocol"""
    same_curve = (P.curve() == Q.curve())
    valid_basis = validate_torsion_basis(P, Q, n)
    # if n is even: both P and Q cannot lie above the (0, 0) point
    not_over_0 = n % 2 == 1 or ((P * (n//2)).xy() != (0, 0) and (Q * (n//2)).xy() != (0, 0))
    return same_curve and valid_basis and not_over_0

def test_tersidh_gen_pub_params():
    print("---: test_tersidh_gen_pub_params()")

    t = 2
    print("t:", t)
    p, A, B, f, _, _ = TerSIDH.gen_pub_params(t)
    print("p:", p)
    print("A:", A)
    print("B:", B)
    print("f:", f)
    assert A == 4 * 5
    assert B == 3 * 7
    assert f == 1
    assert p == A * B * f - 1

    t = 10
    print("t:", t)
    p, A, B, f, _, _ = TerSIDH.gen_pub_params(t)
    print("p:", p)
    print("A:", A)
    print("B:", B)
    print("f:", f)

    t = 50
    print("t:", t)
    p, A, B, f, _, _ = TerSIDH.gen_pub_params(t)
    print("p:", p)
    print("A:", A)
    print("B:", B)
    print("f:", f)

    t = 85
    print("t:", t)
    p, A, B, f, _, _ = TerSIDH.gen_pub_params(t)
    print("p:", p)
    print("A:", A)
    print("B:", B)
    print("f:", f)

class TestT2:

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
        cls.t = 2
        cls.p, cls.A, cls.B, cls.f, _, _ = TerSIDH.gen_pub_params(cls.t)
        cls.n = cls.p + 1

        assert cls.p == 419
        assert cls.A == 20
        assert cls.B == 21
        assert cls.f == 1

        cls.F = GF(cls.p ** 2, modulus=[1,0,1], names="i")
        # Unpack the tuple using comma after i -> ',', same as [0] syntax
        cls.i, = cls.F._first_ngens(1)
        cls.a = 6

        # Montgomery Starting Curve E: y^2 = x^3 + 6x^2 + x
        cls.E = EllipticCurve(cls.F, [0, cls.a, 0, 1, 0])

        if not cls.E.is_supersingular():
            raise ValueError("E is not a supersingular curve")

    @classmethod
    def test_tersidh_generate_kernel_points(cls):
        print("---: test_tersidh_generate_kernel_points()")
        t, i, E, a, n, A = cls.t, cls.i, cls.E, cls.a, cls.n, cls.A

        P = E(124*i + 4, 90*i + 145)
        Q = E(234*i + 290, 175*i + 140)
        PQd = E(75*i + 345, 392*i + 241)

        print(f"xP: {P.x()}")
        print(f"xQ: {Q.x()}")
        print(f"xPQd: {PQd.x()}")
        assert tersidh_base_valid(P, Q, n)
        assert P - Q == PQd

        alice = TerSIDH(t, a, P, Q)
        PA, QA = alice.PA, alice.QA
        print(f"xPA: {PA.x()}")
        print(f"xQA: {QA.x()}")
        assert PA == E(153*i + 357, 375*i + 325)
        assert QA == E(115*i + 221, 414*i + 62)
        assert tersidh_base_valid(PA, QA, A)

        # Generate few "edge-case" secrets

        # 1. s: 0 = '00' -> only KQ gets reduced
        alice = TerSIDH(t, a, P, Q, secret=[0] * t) 
        KP, KQ = alice.kernel_points
        assert KP == alice.PA and KP.order() == 20
        assert KQ == E(0)     and KQ.order() == 1
        print("s: 0")
        print(f"xKP: {KP.x()}")
        print(f"ZKQ: {KQ[2]}")

        # 2. s: 4 = '11' -> only KP gets reduced
        alice = TerSIDH(t, a, P, Q, secret=[1] * t) 
        KP, KQ = alice.kernel_points
        assert KP == E(0)     and KP.order() == 1
        assert KQ == alice.QA and KQ.order() == 20
        print("s: 4")
        print(f"ZKP: {KP[2]}")
        print(f"xKQ: {KQ.x()}")

        # 3. s: 8 = '22' -> both get reduced to neutral element
        alice = TerSIDH(t, a, P, Q, secret=[2] * t)
        KP, KQ = alice.kernel_points
        assert KP == E(0) == KQ
        print("s: 8")
        print(f"ZKP: {KP[2]}")
        print(f"ZKQ: {KQ[2]}")

        # 4. s: 3 = '10' -> KP = [4]PA, KQ = [5]QA 
        alice = TerSIDH(t, a, P, Q, secret=3)
        KP, KQ = alice.kernel_points
        assert KP == alice.PA * 5 and KP.order() == 4
        assert KQ == alice.QA * 4 and KQ.order() == 5
        print("s: 3")
        print(f"xKP: {KP.x()}")
        print(f"xKQ: {KQ.x()}")

        # 5. s: 5 = '12' 
        alice = TerSIDH(t, a, P, Q, secret=5)
        KP, KQ = alice.kernel_points
        print("s: 5")
        assert KP == E(0)         and KP.order() == 1
        assert KQ == alice.QA * 4 and KQ.order() == 5
        print(f"ZKP: {KP[2]}")
        print(f"xKQ: {KQ.x()}")

class TestT15:

    @classmethod
    def setup_class(cls):
        """Prepare variables for running tests"""
        # i is used for constructing fp2 elements: a + b * i

        # p = 419, A = 20, B = 21, f = 1
        cls.t = 15
        cls.p, cls.A, cls.B, cls.f, _, _ = TerSIDH.gen_pub_params(cls.t)
        cls.n = cls.p + 1

        cls.F = GF(cls.p ** 2, modulus=[1,0,1], names="i")
        # Unpack the tuple using comma after i -> ',', same as [0] syntax
        cls.i, = cls.F._first_ngens(1)
        cls.a = 6

        # Montgomery Starting Curve E: y^2 = x^3 + 6x^2 + x
        cls.E = EllipticCurve(cls.F, [0, cls.a, 0, 1, 0])

        if not cls.E.is_supersingular():
            raise ValueError("E is not a supersingular curve")

    @classmethod
    def test_tersidh_state_prepare(cls):
        print("---: test_tersidh_state_prepare()")
        t, i, E, A, n, a = cls.t, cls.i, cls.E, cls.A, cls.n, cls.a
        P = E(29584706441725156430045723882658927461711344409*i + 8723605936412145621728287379910510577218745101, 47578192561637187860896903670848917866094258280*i + 46015566558102970756575555272908862963952112344)
        Q = E(37774940619389439252585951117910616415029401650*i + 49431932972083700521783639676599147784719313148, 30782136400431651285870616173490765995304411624*i + 27192958658720215646695952980861287541095921973)
        PQd = E(27017898316266456978360989490543647986396295282*i + 24459173475706578813044034099488276350833204471, 42936522205384095160941678240474313518339013426*i + 35901259934830632656897115890125102138040843174)

        print(f"xP: {P.x()}")
        print(f"xQ: {Q.x()}")
        print(f"xPQd: {PQd.x()}")
        assert tersidh_base_valid(P, Q, n)

        # secret(ternary): 
        secret = 12722590
        print(f"secret: {secret}")
        alice = TerSIDH(t=t, a=a, P=P, Q=Q, secret=secret, is_bob=False, precalc_public_key=True)
        KP, KQ = alice.kernel_points

        # ter[-1] = 1, means that KQ has order divided by 4, but KP has not
        assert alice.int2ter(secret, t) == [2, 1, 2, 2, 2, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1]
        assert KP.order() == 5 * 11 * 17 * 31 * 47 == 1362295
        assert KQ.order() == 4 * 23 * 41 * 59 * 67 * 103 == 1535803748 
        print(f"ord(KP): {KP.order()}")
        print(f"ord(KQ): {KQ.order()}")

        assert KP == E(34604575071216675628356709894950567972412906170*i + 30198612776533669635778272617725587343553902136, 46845632047117772209135917395701539519507334434*i + 22275638877513543484727541576313459807165920411)
        assert KQ == E(24204599477028428533016801584530225569829265269*i + 8072214545902981320292015944660866527310011603, 59113828676627331816088424668393519835030646708*i + 555295897827375945548498576591575488616350349) 

        # precalc_public_key = True => phi and public_key is already calculated
        assert alice.phi.degree() == KP.order() * KQ.order() == 1362295 * 1535803748
        assert A % alice.phi.degree() == 0

        EA, PB, QB = alice.public_key
        a = EA.a_invariants()[1]

        print(f"aφ(E): {a}")
        print(f"xφ(PB): {PB.x()}")
        print(f"xφ(QB): {QB.x()}")
        assert a == 20823762012321546008743277796759767280614403163*i + 6374112813778050444004853308318914927205109921 
        assert PB == EA(20876846029638301559648453625045466417726680845*i + 6210091974667176375575623846556352273019744322, 44427626681537796372996781492660221933686898778*i + 16957919379089275063011555444131131883737847532)
        assert QB == EA(50225951711665673265080604151407102500355509060*i + 25666354126923489688218829391982119080587425206, 7998706219064803663412177847762019915788171974*i + 40107962334849786454183412795793732388690727370)

    @classmethod
    def test_tersidh_key_exchange(cls):
        print("---: test_tersidh_key_exchange()")
        t, i, E, A, n, a = cls.t, cls.i, cls.E, cls.A, cls.n, cls.a
        # P = E(34165639249113892939227745901234370839942638295*i + 33311322161116326353659629212271304174489381850, 36944833582550965239460623178160990039794984329*i + 36605295699533979934408888983337407290362014912)
        # Q = E(54133067330855987481900044790008400266791298413*i + 26928274597052325448251662784480463387358492758, 31428302797429427800833705878385800309991414470*i + 21252279285061641872803874027595144139655900253)
        # PQd = P - Q

        P   = E(45255132863296035939428643087923170526055812335*i + 35207532789640029607392085315164843785886696913, 62549254984344532813993825360042331319918815620*i + 53588053495568172647974770086454463944747703409) 
        Q   = E(62188135383560125911606431706677411561756802948*i + 52478616152221885238374224345805957897724858098, 34710578818202657047452235142231297688646060235*i + 52458580608776717800219352585205182498132453796)
        PQd = E(31403620116220219651357966569215397638854000763*i + 20465179760444544011039140556083357241775723149, 40807487265125135025757204537726279318148430294*i + 29413969703366854039155026095354974735052008388)
        print(f"xP: {P.x()}")
        print(f"xQ: {Q.x()}")
        print(f"xPQd: {PQd.x()}")
        assert PQd == P - Q
        assert tersidh_base_valid(P, Q, n)

        a_secret = 6631513
        print(f"a_secret: {a_secret}")
        alice = TerSIDH(t, a, P=P, Q=Q, secret=a_secret, is_bob=False, precalc_public_key=True)

        b_secret = 4980130 
        print(f"b_secret: {b_secret}")
        bob   = TerSIDH(t, a, P=P, Q=Q, secret=b_secret,  is_bob=True, precalc_public_key=True)

        EA, PBA, QBA = alice.public_key
        aEA = EA.a_invariants()[1]
        print(f"aφA(E0): {aEA}")
        print(f"xφA(PB): {PBA.x()}")
        print(f"xφA(QB): {QBA.x()}")
        assert aEA == 54801376089903970925339307196694608662317200341*i + 3248665947417223701882116168791158651436364321
        assert PBA == EA(43868479477697879809566585639320978671149440254*i + 7059168512921389348457558705746610439407082700, 40033696492354779234469600977932461354581402526*i + 41484416941716509897288416565024722787722837763)
        assert QBA == EA(59025181202054368531235216709432152582925746866*i + 32671210671678364256011079742850149244716416643, 30908442585584827037329478060849336947348599785*i + 33299724605681349170912816283549452722086975612)

        EB, PAB, QAB = bob.public_key
        aEB = EB.a_invariants()[1]
        print(f"aφB(E0): {aEB}")
        print(f"xφB(PA): {PAB.x()}")
        print(f"xφB(QA): {QAB.x()}")
        assert aEB == 3181532679381964117351769805003588383567570283*i + 8002906843868878916100601361840467966138396994
        assert PAB == EB(19693539544337260969862178248809154682696090848*i + 32057768515926752112257154360884331701455434315, 8706402505749026696092462919815217550192784942*i + 52096932948786443913508572791963895616147378893)
        assert QAB == EB(49661089078839342568423577829915630875462416634*i + 52781414101202328567169277344100576338432492581, 36972311997144233498723792265423430315509744335*i + 29885959950124116988283923902516211043279876169)

        j_inv1 = alice.key_exchange(EB, PAB, QAB)
        j_inv2 = bob.key_exchange(EA, PBA, QBA)
        print(f"j(EBA): {j_inv1}")
        print(f"j(EAB): {j_inv2}")
        assert j_inv1 == j_inv2 == 16477822473601326380854754948643703255616912674*i + 12600493404726659034043219507047767318775160385

def main():
    test_tersidh_gen_pub_params()

    # Small TerSIDH t = 2
    TestT2.setup_class()
    TestT2.test_tersidh_generate_kernel_points()

    TestT15.setup_class()
    TestT15.test_tersidh_state_prepare()
    TestT15.test_tersidh_key_exchange()

if __name__ == "__main__":
    main()





