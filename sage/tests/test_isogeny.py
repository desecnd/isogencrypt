#!/usr/bin/sage -python

from sage.all import GF, EllipticCurve
from isogencrypt_sage.isogeny import sample_torsion_basis_smooth

class TestIsogeny:

    def setup_method(self):
        # 2^3 * 103 - 1
        self.p = 823

        self.F = GF(self.p**2, names=('i',), modulus=[1,0,1])
        (self.i,) = self.F._first_ngens(1)
        self.E = EllipticCurve(self.F, [0, 6, 0, 1, 0])
        assert self.E.is_supersingular()

    def test_random_sampling(self):
        r = self.p + 1

        def above_zero(R) -> bool:
            return (R * (r//2)).x() == 0

        for _ in range(10):

            P, Q = sample_torsion_basis_smooth(self.E, r, 'P')
            assert above_zero(P) and not above_zero(Q)
            # assert (P * (r//2)).x() == 0 and (Q * (r//2)).x() != 0

            P, Q = sample_torsion_basis_smooth(self.E, r, 'Q')
            assert not above_zero(P) and above_zero(Q)
            # assert (Q * (r//2)).x() == 0 and (P * (r//2)).x() != 0

            P, Q = sample_torsion_basis_smooth(self.E, r, 'none')
            assert not above_zero(P) and not above_zero(Q)
            # assert (P * (r//2)).x() != 0 and (Q * (r//2)).x() != 0


