#/usr/bin/sage 

from sage.all import EllipticCurve, Primes, randint, order_from_multiple, gcd, is_prime, prod, GF
from verifiers.isogeny import sample_quadratic_root_of_unity, sample_torsion_basis_smooth, mont_coef, check_points_torsion_basis

class MSIDH:
    def __init__(self, p: int, A: int, B: int, E0, P = None, Q = None, secret: int = None, mask: int = None, is_bob: bool = False, mont_model: bool = False):
        self.p = p
        self.A = A
        self.B = B
        self.is_bob = is_bob

        self.E0 = E0
        self.P = P
        self.Q = Q
        self.secret = secret
        self.mask = mask
        self.mont_model = mont_model

        # Characteristic of the starting Curve must match
        assert E0.base().characteristic() == p

        if is_bob:
            self.name = "Bob"
            self.A, self.B = self.B, self.A
        else:
            self.name = "Alice"

        if self.P is None or self.Q is None:
            print(f"{self.name}: Sampling Torsion Basis (P, Q)...")
            self.P, self.Q = sample_torsion_basis_smooth(self.E0, self.A)
        else:
            assert self.P.curve() == self.E0
            assert self.Q.curve() == self.E0

        if self.mask is None:
            print(f"{self.name}: Sampling mask...")
            self.mask = sample_quadratic_root_of_unity(self.B)
        else:
            # assert pow(self.mask, 2, self.B) == 1
            pass
        
        if self.secret is None:
            print(f"{self.name}: Generating secret...")
            self.secret = randint(0, self.A)
        else:
            assert self.secret in range(self.A)

        self.phi_ker = self.P + self.secret * self.Q
        if self.mont_model:
            A_ = mont_coef(self.phi_ker)
            E_ = EllipticCurve(E0.base(), [0, A_, 0, 1, 0])
            self.phi = self.E0.isogeny(self.phi_ker, codomain=E_, algorithm='factored')
        else:
            self.phi = self.E0.isogeny(self.phi_ker, algorithm='factored')
        self.E = self.phi.codomain()
        self.key = None
    
    @property
    def basis(self):
        return (self.P, self.Q)

    @classmethod
    def gen_pub_params(cls, t: int):
        # Number of primes to multiply  ~ t/2 bit classical and ~t/4 quantum security
        P = Primes()
        ll = [ P.unrank(2 * i) for i in range((t+1)//2) ]
        qq = [ P.unrank(2 * i + 1) for i in range(t//2) ]
        ll[0] = 4

        assert len(ll) + len(qq) == t
        A = prod(ll)
        B = prod(qq)
        AB = A * B
        assert gcd(A, B) == 1

        for f in range(1000):
            p = f * AB - 1
            if is_prime(p):
                break
        else:
            raise ValueError("Cannot find cofactor for p")

        return p, A, B, f

    def gen_pubkey(self, PB, QB) -> tuple:
        assert PB.curve() == self.E0
        assert QB.curve() == self.E0
        assert check_points_torsion_basis(PB, QB, self.B)

        # Apply masking
        PB = self.mask * self.phi(PB)
        QB = self.mask * self.phi(QB)
        return (self.E, PB, QB)

    def key_exchange(self, EB, BPA, BQA):
        assert BPA.curve() == EB
        assert BQA.curve() == EB 
        # Weil pairing property of the curve
        assert BPA.weil_pairing(BQA, self.A) == self.P.weil_pairing(self.Q, self.A) ** self.B

        self.tau_ker = BPA + BQA * self.secret 
        if self.mont_model:
            A_ = mont_coef(self.tau_ker)
            E_ = EllipticCurve(self.E0.base(), [0, A_, 0, 1, 0])
            self.tau = EB.isogeny(self.tau_ker, codomain=E_, algorithm="factored")
        else:
            self.tau = EB.isogeny(self.tau_ker, algorithm="factored")

        self.EK = self.tau.codomain()
        return self.EK.j_invariant()

if __name__ == '__main__':
    p, A, B, f = MSIDH.gen_pub_params(6)

    F = GF(p**2, names=('i',))
    (i,) = F._first_ngens(1)

    E0 = EllipticCurve(F, [0, 6, 0, 1, 0])

    Alice = MSIDH(p, A, B, E0, is_bob=False)
    Bob = MSIDH(p, A, B, E0, is_bob=True)

    pubkey_alice = Alice.gen_pubkey(*Bob.basis)
    pubkey_bob = Bob.gen_pubkey(*Alice.basis)

    jinv_alice = Alice.key_exchange(*pubkey_bob)
    jinv_bob = Bob.key_exchange(*pubkey_alice)
    print(jinv_alice == jinv_bob)