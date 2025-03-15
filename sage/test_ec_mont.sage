p = 431
F.<i> = GF(p^2, modulus=[1,0,1])
# Montgomery Starting Curve
# E: y^2 = x^3 + 6x^2 + x
E = EllipticCurve(F, [0, 6, 0, 1, 0])


def test_formula_xADD():
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