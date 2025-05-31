#!/usr/bin/sage -python

import pytest

from sage.all import GF, EllipticCurve
from isogencrypt_sage.isogeny import sample_torsion_basis_smooth, mont_coef, mont_isog

class TestIsogeny:

    def test_random_sampling(self):
        # 2^3 * 103 - 1
        p = 823
        F = GF(p**2, names=('i',), modulus=[1,0,1])
        (i,) = F._first_ngens(1)
        E = EllipticCurve(F, [0, 6, 0, 1, 0])
        assert E.is_supersingular()
        r = p + 1

        def above_zero(R) -> bool:
            return (R * (r//2)).x() == 0

        for _ in range(10):
            P, Q = sample_torsion_basis_smooth(E, r, 'P')
            assert above_zero(P) and not above_zero(Q)

            P, Q = sample_torsion_basis_smooth(E, r, 'Q')
            assert not above_zero(P) and above_zero(Q)

            P, Q = sample_torsion_basis_smooth(E, r, 'none')
            assert not above_zero(P) and not above_zero(Q)

    def test_ISOG2_bad_point_error(self):
        # p = 4 * 5 * 7 - 1
        p = 139
        F = GF(p ** 2, modulus=[1,0,1], names="i")
        i, = F._first_ngens(1)

        # Montgomery Starting Curve E: y^2 = x^3 + 6x^2 + x
        E = EllipticCurve(F, [0, 6, 0, 1, 0])
        assert E.is_supersingular()

        K2 = E(0, 0)
        assert K2.order() == 2
        print(f"xK: {K2.x()}")

        # Formula for obtaining the A' for 2-isogeny
        # This will give value 2 which is incorrect (singular curve)
        A2 = 2 * (1 - 2 * K2.x() ** 2)

        # _mont_coef_2 checks for invalid kernel
        with pytest.raises(ValueError):
            mont_coef(K2)

        # Montgomery curve is singular for A^2 = 4
        assert A2**2 == 4

        # This must fail due to A2 being equal to +-2
        # Curve is Singular
        with pytest.raises(ArithmeticError):
            EllipticCurve(F, [0, A2, 0, 1, 0])

