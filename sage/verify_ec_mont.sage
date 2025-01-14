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

def verify_test_small_xADD():

    P = E(271*i + 259, 422*i + 97)
    Q = E(335*i + 262, 69*i + 198)
    PQdiff = E(411*i + 143, 245*i + 213)
    assert PQdiff == P - Q
    PQsum = P + Q

    assert PQsum == E(106*i + 416 , 111*i + 405)
    assert PQsum.x() == 416 + 106 * i

def verify_test_small_xLADDER3PT():
    P = E(271*i + 259, 422*i + 97)
    Q = E(335*i + 262, 69*i + 198)
    n = 87

    R = P + n * Q
    assert R == E(45*i + 360, 249*i + 429)
    assert R.x() == 360 + 45 * i


if __name__ == '__main__':

    verify_test_small_xDBL()
    verify_test_small_xADD()
    verify_test_small_xLADDER3PT()


