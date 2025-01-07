# p = 2^4 * 3^3 - 1
p = 431
F.<i> = GF(p**2, modulus=[1,0,1])

# Montgomery Starting Curve
# E: y^2 = x^3 + 6x^2 + x
E = EllipticCurve(F, [0, 6, 0, 1, 0])

# Point with order 432
P = E((292 + 15 * i, 281 + 235 * i))
# Q = [2]P
Q = E((61 + 184 * i, 395 + 90 * i))
assert Q == 2 * P

