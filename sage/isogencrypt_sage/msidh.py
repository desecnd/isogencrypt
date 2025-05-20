#/usr/bin/sage -python

from dataclasses import dataclass, asdict
import json

from sage.all import EllipticCurve, Primes, randint, gcd, is_prime, prod, GF
from isogencrypt_sage.isogeny import sample_quadratic_root_of_unity, sample_torsion_basis_smooth, mont_coef, validate_torsion_basis

class MSIDH:
    def __init__(self, p: int, A: int, B: int, E0, P = None, Q = None, secret: int | None = None, mask: int | None = None, is_bob: bool = False, mont_model: bool = False):
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
            # print(f"{self.name}: Sampling Torsion Basis (P, Q)...")
            self.P, self.Q = sample_torsion_basis_smooth(self.E0, self.A)
        else:
            assert self.P.curve() == self.E0
            assert self.Q.curve() == self.E0
            assert validate_torsion_basis(self.P, self.Q, self.A)

        if self.mask is None:
            # print(f"{self.name}: Sampling mask...")
            self.mask = sample_quadratic_root_of_unity(self.B)
        else:
            # assert pow(self.mask, 2, self.B) == 1
            pass
        
        if self.secret is None:
            # print(f"{self.name}: Generating secret...")
            self.secret = randint(0, self.A)
        else:
            assert self.secret in range(self.A)

        self.phi_ker = self.P + self.secret * self.Q
        self.phi_ker.set_order(self.A)
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
        P = Primes(proof=False)
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
        assert validate_torsion_basis(PB, QB, self.B)

        # Apply masking
        PB = self.mask * self.phi(PB)
        QB = self.mask * self.phi(QB)
        return (self.E, PB, QB)

    def key_exchange(self, EB, BPA, BQA):
        assert BPA.curve() == EB
        assert BQA.curve() == EB 
        # Weil pairing property of the curve
        assert BPA.weil_pairing(BQA, self.A) == self.P.weil_pairing(self.Q, self.A) ** self.B # type: ignore

        self.tau_ker = BPA + BQA * self.secret 
        self.tau_ker.set_order(self.A)
        if self.mont_model:
            A_ = mont_coef(self.tau_ker)
            E_ = EllipticCurve(self.E0.base(), [0, A_, 0, 1, 0])
            self.tau = EB.isogeny(self.tau_ker, codomain=E_, algorithm="factored")
        else:
            self.tau = EB.isogeny(self.tau_ker, algorithm="factored")

        self.EK = self.tau.codomain()
        return self.EK.j_invariant()

@dataclass
class MSIDHBenchTask:
    t: int
    a: str
    xP: str
    yP: str
    xQ: str
    yQ: str
    xR: str
    yR: str
    
    @classmethod
    def from_dict(cls, obj: dict):
        return MSIDHBenchTask(**obj)

    def to_dict(self) -> dict:
        return asdict(self) 

def load_msidh_bench_tasks(filename: str) -> list[MSIDHBenchTask]:
    with open(filename) as in_file:
        msidh_objs = json.loads(in_file.read())
    
    if "bench_tasks" not in msidh_objs.keys():
        raise ValueError("Loaded JSON object from MSIDH bench tasks does not contain 'bench_tasks' key")

    msidh_objs = msidh_objs["bench_tasks"]
    msidh_bts = [ MSIDHBenchTask.from_dict(obj) for obj in msidh_objs ]
    return msidh_bts

def store_msidh_bench_tasks(filename: str, msidh_bts: list[MSIDHBenchTask]): 
    msidh_objs = { 
        "bench_tasks": [bt.to_dict() for bt in msidh_bts ]
    }
    with open(filename, 'w') as out_file:
        out_file.write(json.dumps(msidh_objs, indent=4))

if __name__ == '__main__':
    p, A, B, f = MSIDH.gen_pub_params(20)

    F = GF(p**2, names=('i',), modulus=[1, 0, 1])
    (i,) = F._first_ngens(1)

    E0 = EllipticCurve(F, [0, 6, 0, 1, 0])
    assert E0.is_supersingular()

    Alice = MSIDH(p, A, B, E0, is_bob=False, mont_model=True)
    Bob = MSIDH(p, A, B, E0, is_bob=True, mont_model=True)

    pubkey_alice = Alice.gen_pubkey(*Bob.basis)
    pubkey_bob = Bob.gen_pubkey(*Alice.basis)

    jinv_alice = Alice.key_exchange(*pubkey_bob)
    jinv_bob = Bob.key_exchange(*pubkey_alice)
    assert jinv_alice == jinv_bob