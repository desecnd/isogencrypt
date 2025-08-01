#/usr/bin/sage
from sage.all import EllipticCurve, factor, order_from_multiple, Primes, prod, gcd, is_prime, GF,  randint, set_random_seed
from isogencrypt_sage.isogeny import sample_torsion_basis_smooth, validate_torsion_basis, mont_isog
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_finite_field as Point 

class TerSIDH:

    MAX_DIGIT = 2

    @classmethod
    def gen_pub_params(cls, t: int) -> tuple[int, int, int, int, list[int], list[int]]:
        """Constructs public parameters for TerSIDH protocol

        Args:
            t: Parameter related to security (in bits) - number of prime numbers per side.

        Raises:
            ValueError: Cannot find prime number in 1000 iterations
        """

        P = Primes(proof=False)
        A_primes = [ P.unrank(2 * i) for i in range(t) ]
        B_primes = [ P.unrank(2 * i + 1) for i in range(t) ]
        # Use 4 instead of 2 to always generate prime = 3 mod 4
        A_primes[0] = 4

        assert len(A_primes) + len(B_primes) == 2 * t

        A_deg = prod(A_primes)
        B_deg = prod(B_primes)
        AB = A_deg * B_deg
        assert gcd(A_deg, B_deg) == 1

        for f in range(1000):
            p = f * AB - 1
            if is_prime(p):
                break
        else:
            raise ValueError("Cannot find cofactor for p")

        return p, A_deg, B_deg, f, A_primes, B_primes

    @classmethod
    def generate_kernel_coeffs(cls, secret: list[int], primes: list[int]) -> tuple[int, int, int, int]:
        cP, nP = 1, 1
        cQ, nQ = 1, 1
        for digit, p in zip(secret, primes):
            if digit == 0:
                nP *= p
                cQ *= p
            elif digit == 1:
                nQ *= p
                cP *= p
            elif digit == 2:
                # Decrease the order of the isogeny
                cP *= p
                cQ *= p
            else:
                raise ValueError("Wrong Digit Value")

        return cP, nP, cQ, nQ

    @property
    def is_alice(self) -> bool:
        return not self.is_bob

    def __init__(self, t: int, a, P = None, Q = None, secret: list[int] | str | None = None, mask: int | None = None, is_bob: bool = False, precalc_public_key: bool = False):
        """Initialize TerSIDH interface. All variables are defines from Alice perspective."""

        p, A_deg, B_deg, f, A_primes, B_primes = self.gen_pub_params(t)

        A_deg, B_deg = (B_deg, A_deg) if is_bob else (A_deg, B_deg)
        A_primes, B_primes = (B_primes, A_primes) if is_bob else (A_primes, B_primes)

        F = GF(p**2, names=('i',), modulus=[1, 0, 1])
        (i,) = F._first_ngens(1)
        n = p + 1 

        # Montgomery Curve Model with coeff a near x^2
        E0 = EllipticCurve(F, [0, a, 0, 1, 0])
        assert E0.order() == n**2
        assert E0.is_supersingular()

        if P is None or Q is None:
            # Alice cannot use torsion basis, which is not 'none' category - even isogeny will not work
            P, Q = sample_torsion_basis_smooth(E0, n, 'none')
        else:
            P = E0(P)
            Q = E0(Q)

        assert validate_torsion_basis(P, Q, n)

        PA = P * (n // A_deg)
        QA = Q * (n // A_deg)

        PB = P * (n // B_deg)
        QB = Q * (n // B_deg)

        # if we are Alice, none of the points can lie above point E(0, 0)
        if not is_bob:
            assert PA * (A_deg // 2) != E0(0, 0)
            assert QA * (A_deg // 2) != E0(0, 0)

        # Generate point Masking information
        if mask is None:
            while mask := randint(1, B_deg):
                if gcd(mask, B_deg) == 1:
                    break

        # Generate Ternary Secret
        if secret is None:
            secret = [ randint(0, self.MAX_DIGIT) for _ in range(t) ] 
        else:
            if isinstance(secret, str):
                secret = [ int(x) for x in secret ]
            elif not isinstance(secret, list):
                raise ValueError("Incorrect secret type")

            assert all(x in range(0, self.MAX_DIGIT + 1) for x in secret)
            assert len(secret) == t

        self.cP, self.nP, self.cQ, self.nQ = self.generate_kernel_coeffs(secret, A_primes)
        self.is_bob = is_bob
        self.secret = secret
        self.mask = mask
        self.A_deg = A_deg
        self.P, self.Q = P, Q
        self.E0 = E0
        self.PA, self.QA = PA, QA
        self.PB, self.QB = PB, QB
        self.E_pub, self.P_pub, self.Q_pub = None, None, None

        if precalc_public_key:
            self.prepare_public_key()

    @property
    def start_basis(self) -> tuple[object, Point, Point]:
        return (self.E0, self.P, self.Q)

    @property
    def public_key(self) -> tuple[object, Point, Point]:
        if self.E_pub is None or self.P_pub is None or self.Q_pub is None:
            raise ValueError("Public basis was not calculated yet. Call self.prepare() first.")

        return (self.E_pub, self.P_pub, self.Q_pub)

    def prepare_public_key(self) -> tuple:

        KP = self.PA * self.cP
        assert KP.order() == self.nP
        KQ = self.QA * self.cQ
        assert KQ.order() == self.nQ

        phi_P = mont_isog(KP)
        # Send the point through the isogeny
        KQ = phi_P(KQ)

        phi_Q = mont_isog(KQ)
        EQ = phi_Q.codomain()

        self.E_pub = EQ
        self.P_pub = self.mask * phi_Q(phi_P(self.PB))
        self.Q_pub = self.mask * phi_Q(phi_P(self.QB))
        # precalculated above
        return self.public_key


    def key_exchange(self, EB_pub, PB_pub: Point, QB_pub: Point):
        # Verify that we have obtained our new torsion basis that is valid
        assert validate_torsion_basis(PB_pub, QB_pub, self.A_deg)
        assert PB_pub.curve() == EB_pub
        assert QB_pub.curve() == EB_pub

        KP = PB_pub * self.cP
        assert KP.order() == self.nP
        KQ = QB_pub * self.cQ
        assert KQ.order() == self.nQ

        phi_P = mont_isog(KP)
        KQ = phi_P(KQ)
        phi_Q = mont_isog(KQ)
        return phi_Q.codomain().j_invariant()
        
if __name__ == "__main__":
    set_random_seed(0)
    t = 30
    a = 6
    print("Init Alice")
    alice = TerSIDH(t, a)
    _, P, Q = alice.start_basis
    print("Init Bob")
    bob = TerSIDH(t, a, P, Q, is_bob=True)

    # Generate public keys for Alice and Bob
    print("Prepare Alice")
    alice.prepare_public_key()
    print("Prepare Bob")
    bob.prepare_public_key()

    print("Exchange Alice")
    ss_a = alice.key_exchange(*bob.public_key)
    print("Exchange Bob")
    ss_b = bob.key_exchange(*alice.public_key)
    assert ss_a == ss_b
    print("Ok")

