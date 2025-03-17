#!/usr/bin/sage

def find_ell_prime_number(factors: list[int]) -> int:
    """Find prime number that contains specific factors under of p + 1 and makes a valid supersingular ec graph"""

    ITERS = 100
    for i in range(1, ITERS):
        p = prod(factors) * 2 * i - 1
        if not is_prime(p):
            continue
        try:
            # Modulus must be x^2 + 1, otherwise the arithmetic may be wrong
            U.<u> = GF(p^2, modulus=[1,0,1])
            if (EllipticCurve(U, [1, 0]).is_supersingular()):
                return p
        except:
            continue
    else:
        raise ValueError(f"Cannot find valid prime number in {ITERS} iterations.")

def test_formula_xADD():
    """Test that the formula used for xADD works correctly"""

    p = 431
    F.<i> = GF(p^2, modulus=[1,0,1])
    # Montgomery starting curve E: y^2 = x^3 + 6x^2 + x
    E = EllipticCurve(F, [0, 6, 0, 1, 0])

    print("TEST ---: xADD")
    P = E(107*i + 11, 301*i + 428)
    Q = E(71*i + 243, 61*i + 18)

    R = P - Q 
    assert R == E(286*i + 119, 301*i + 124)

    T = P + Q
    assert T == E(62*i + 260, 42*i + 180)

    print(f"x(P + Q) = {T.x()}")

    print("1) (X : Z) = (x : 1)")
    XP = P.x()
    XQ = Q.x() 
    XR = R.x()
    ZP = ZQ = ZR = 1

    # xADD formula in projective coordinates: P + Q = (X' : Z')
    # X' = ZR * (XP * XQ - ZP * ZQ)^2 
    # Z' = XR * (XP * ZQ - XQ * ZP)^2
    X = ZR * (XP * XQ - ZP * ZQ)^2
    Z = XR * (XP * ZQ - XQ * ZP)^2

    print(f"{X = }")
    print(f"{Z = }")
    x = X/Z
    print(f"> X/Z = {x = }")
    assert x == T.x()

    print("2) (X : Z) = (ax : a)")
    # Test in projective space
    # Multiply by 'random' scalars => should still work
    XP, ZP = [ (416*i + 7)   * coord for coord in [XP, ZP] ]
    XQ, ZQ = [ (247*i + 159) * coord for coord in [XQ, ZQ] ]
    XR, ZR = [ (75*i + 79)   * coord for coord in [XR, ZR] ]

    X = ZR * (XP * XQ - ZP * ZQ)^2
    Z = XR * (XP * ZQ - XQ * ZP)^2
    print(f"{X = }")
    print(f"{Z = }")
    x = X/Z
    print(f"> X/Z = {x = }")
    assert x == T.x()


if __name__ == '__main__':
    test_formula_xADD()