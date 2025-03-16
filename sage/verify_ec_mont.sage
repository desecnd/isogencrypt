#!/usr/bin/sage 

p = 431
F.<i> = GF(p^2, modulus=[1,0,1])
# Montgomery Starting Curve
# E: y^2 = x^3 + 6x^2 + x
E = EllipticCurve(F, [0, 6, 0, 1, 0])

def verify_test_small_xDBL():
    # Point with order 432
    P = E((292 + 15 * i, 281 + 235 * i))
    Q = E(61 + 184 * i, 395 + 90 * i)

    # Q = [2]P
    assert 2 * P == Q
    assert Q.x() == 61 + 184 * i

# Function for testing xADD:
def verify_test_xADD_small():
    print("---: test_xADD_small")

    P = E(271*i + 259, 422*i + 97)
    print(f"xP: {P.x()}")

    Q = E(335*i + 262, 69*i + 198)
    print(f"xQ: {Q.x()}")

    PQdiff = E(411*i + 143, 245*i + 213)
    assert PQdiff == P - Q

    print(f"xP-Q: {PQdiff.x()}")
    PQsum = P + Q

    assert PQsum == E(106*i + 416 , 111*i + 405)
    print(f"xP+Q: {PQsum.x()}")

    assert PQsum.x() == 416 + 106 * i

def verify_test_criss_cross_small():
    print("---: test_criss_cross_small")
    x = F(416*i + 175)
    y = F(112*i + 179)
    z = F(235*i + 107)
    w = F(183*i + 197)

    e = x * w + y * z
    f = x * w - y * z
    # assert e == 206*i + 256
    # assert f == 180*i + 52 
    print(x * w)
    print(y * z)

    print(f"x: {x}")
    print(f"y: {y}")
    print(f"z: {z}")
    print(f"w: {w}")
    print(f"xw+yz: {e}")
    print(f"xw-yz: {f}")

    
def verify_test_xLADDER3PT_small():
    print("---: test_xLADDER3PT_small")
    P = E(271*i + 259, 422*i + 97)
    print(f"xP: {P.x()}")
    Q = E(335*i + 262, 69*i + 198)
    print(f"xQ: {Q.x()}")
    PQd = P - Q
    print(f"xP-Q: {PQd.x()}")
    n = 87
    print(f"n: {n}")

    R = P + n * Q
    print(f"x(P+nQ): {R.x()}")
    assert R == E(45*i + 360, 249*i + 429)

    T = P - n * Q
    print(f"x(P-nQ): {T.x()}")

def verify_point_normalize_coords():
    print("---: point_normalize_coords()")

    X = 395*i + 201
    print(f"X: {X}")
    Y = 272*i + 286
    print(f"Y: {Y}")

    X_ = X / Y
    print(f"X': {X_}")
    assert X_ == 12*i + 95


if __name__ == '__main__':

    verify_test_small_xDBL()
    verify_test_xADD_small()
    verify_test_criss_cross_small()
    verify_test_xLADDER3PT_small()
    verify_point_normalize_coords()


