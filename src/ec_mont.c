#include <gmp.h>

#include "ec_mont.h"

static int g_layer_ctx_initialized = 0;
static fp2_t lr1, lr2, lr3, lr4;

void ec_layer_ctx_init() {
    if (g_layer_ctx_initialized) 
        return;
    fp2_init(&lr1);
    fp2_init(&lr2);
    fp2_init(&lr3);
    fp2_init(&lr4);
    g_layer_ctx_initialized = 1;
}

void ec_layer_ctx_clear() {
    if (!g_layer_ctx_initialized)
        return;
    fp2_clear(&lr1);
    fp2_clear(&lr2);
    fp2_clear(&lr3);
    fp2_clear(&lr4);
    g_layer_ctx_initialized = 0;
}

void xDBL(point_t R, const point_t P, const fp2_t A24p, const fp2_t C24) {

    // P = (X : Z)
    fp2_sub(lr1, P->X, P->Z); // t0 = X - Z
    fp2_add(lr2, P->X, P->Z); // t1 = X + Z

    // Warning! Must use "safe" function for square-ing
    fp2_sq_safe(lr1); // t0 = (X - Z)^2
    fp2_sq_safe(lr2); // t1 = (X + Z)^2

    fp2_mul_unsafe(R->Z, lr1, C24);  // Z' = (X - Z)^2 * C24
    fp2_mul_unsafe(R->X, R->Z, lr2); // X' = (X + Z)^2 * (X - Z)^2 * C24

    fp2_sub(lr2, lr2, lr1);          // t1 = (X + Z)^2 - (X - Z)^2
    fp2_mul_unsafe(lr1, A24p, lr2); // t0 = A24p * [(X + Z)^2 - (X - Z)^2]

    // Z' = A24p * [(X + Z)^2 - (X - Z)^2] + (X - Z)^2 * C24
    fp2_add(R->Z, R->Z, lr1);
    // Z' = { A24p * [(X + Z)^2 - (X - Z)^2] + C24 (X - Z)^2] } [(X + Z)^2 - (X
    // - Z)^2]
    fp2_mul_safe(R->Z, lr2);

}

void xDBLe(point_t R, const point_t P, const fp2_t A24p, const fp2_t C24,
           const int e) {
    point_set(R, P);
    // Repeat the step of doubling multiple times
    for (int i = 0; i < e; i++) {
        xDBL(R, R, A24p, C24);
    }
}

void xADD(point_t PQsum, const point_t P, const point_t Q,
          const point_t PQdiff) {
    // Function is argument safe for calling PQSum = Q or P or PQdiff
    // Given P, Q and PQdiff = P-Q output P+Q

    // xADD formula in projective coordinates: P + Q = (X' : Z')
    // X' = ZR * (XP * XQ - ZP * ZQ)^2
    // Z' = XR * (XP * ZQ - XQ * ZP)^2
    //
    // Optimised version:
    // X' = ZR * [(XP - ZP)(XQ + ZQ) + (XP + ZP)(XQ - ZQ)]^2
    // Z' = XR * [(XP - ZP)(XQ + ZQ) - (XP + ZP)(XQ - ZQ)]^2

    // Defining additional variables:
    // a: XP + ZP
    // b: XP - ZP
    // c: XQ + ZQ
    // d: XQ - ZQ
    //
    // X' = ZR * (b * c + a * d)^2
    // Z' = XR * (b * c - a * d)^2
    fp2_add(lr1, P->X, P->Z); // t0 = a: XP + ZP
    fp2_sub(lr2, P->X, P->Z); // t1 = b: XP - ZP
    fp2_add(lr3, Q->X, Q->Z); // t2 = c: xQ + zQ
    fp2_sub(lr4, Q->X, Q->Z); // t3 = d: xQ - zQ

    fp2_mul_safe(lr1, lr4); // t0 = t0 * t3: a * d
    fp2_mul_safe(lr2, lr3); // t1 = t1 * t2: b * c

    fp2_add(lr3, lr1, lr2); // t2 = t0 + t1: ad + bc
    fp2_sub(lr4, lr1, lr2); // t3 = t0 - t1: ad - bc

    fp2_sq_unsafe(lr1, lr3); // t0 = t2^2: (ad + bc)^2
    fp2_sq_unsafe(lr2, lr4); // t1 = t3^2: (ad - bc)^2

    // We must use t2 to not override the used later PQdiff->X in scenario when
    // PQsum = PQdiff t2 = t0 * Z(P-Q): (ad + bc)^2 * Z(P-Q)
    fp2_mul_unsafe(lr3, lr1, PQdiff->Z);
    // Z' = t1 * X(P-Q): (ad - bc)^2 * Z(P-Q)
    fp2_mul_unsafe(PQsum->Z, lr2, PQdiff->X);
    // X' = t2
    fp2_set(PQsum->X, lr3);
    // TODO: fp2_swap?
}

// Out: P = 2P, Q = P + Q
// TODO: check argument-safeness
// TODO: think about expressing function
// as set of arithmetic instructions
// instead of combinatin of 2 function calls
inline void xDBLADD(point_t P, point_t Q, const point_t PQdiff, const fp2_t A24p,
             const fp2_t C24) {
    xADD(Q, P, Q, PQdiff);
    xDBL(P, P, A24p, C24);
}

void xLADDER(point_t R0, const point_t P, const mpz_t m, const fp2_t A24p,
             const fp2_t C24) {
    assert(mpz_sgn(m) > 0 && "Given scalar m must be nonnegative");

    mpz_t r;
    mpz_init(r);

    point_t R1;
    point_init(&R1);

    // R0 = P, R1 = [2]R
    point_set(R0, P);
    xDBL(R1, P, A24p, C24);

    // Get number of "active" bits
    int n_bits = mpz_sizeinbase(m, 2);

    // Iterate over bits downwards (leading bit not included).
    // Invariant of the algorithm R1 - R0 = P
    for (int bit = n_bits - 2; bit >= 0; bit--) {
        // Bit is equal to 1
        if (mpz_tstbit(m, bit)) {
            // R1 = [2]R1; R0 = R0 + R1
            xDBLADD(R1, R0, P, A24p, C24);
        } else {
            // R0 = [2]R0; R1 = R0 + R1
            xDBLADD(R0, R1, P, A24p, C24);
        }
    }

    mpz_clear(r);
    point_clear(&R1);
}

void xLADDER_int(point_t R0, const point_t P, long int m, const fp2_t A24p,
                 const fp2_t C24) {
    assert(m > 0 && "Given scalar m must be nonnegative");

    point_t R1;
    point_init(&R1);

    // R0 = P, R1 = [2]R
    point_set(R0, P);
    xDBL(R1, P, A24p, C24);

    // Get number of "active" bits (count until leading bit is found)
    int bits = 0;
    for (long int x = m; x > 0; x /= 2)
        bits++;

    // Iterate over bits downwards (leading bit not included).
    // Invariant of the algorithm R1 - R0 = P
    for (int bit = bits - 2; bit >= 0; bit--) {

        // bit = 1
        if (m & (1 << bit)) {
            // R1 = [2]R1; R0 = R0 + R1
            xDBLADD(R1, R0, P, A24p, C24);
        } else {
            // R0 = [2]R0; R1 = R0 + R1
            xDBLADD(R0, R1, P, A24p, C24);
        }
    }

    point_clear(&R1);
}

// calculate P = P + [m]Q
void xLADDER3PT_int(point_t P, point_t Q, point_t PQdiff, long int m,
                    const fp2_t A24p, const fp2_t C24) {
    assert(m >= 0 && "Given scalar m must be nonnegative");

    while (m > 0) {
        if (m & 1)
            xDBLADD(Q, P, PQdiff, A24p, C24);
        else
            xDBLADD(Q, PQdiff, P, A24p, C24);
        m /= 2;
    }
}

// calculate P = P + [m]Q
void xLADDER3PT(point_t P, point_t Q, point_t PQdiff, const mpz_t m,
                const fp2_t A24p, const fp2_t C24) {
    assert(mpz_sgn(m) >= 0 && "Given scalar m must be nonnegative");

    mpz_t r, n;
    mpz_init(r);
    mpz_init_set(n, m);

    while (mpz_sgn(n) > 0) {
        // m = m//2; r = m % 2
        mpz_fdiv_qr_ui(n, r, n, 2);

        if (mpz_sgn(r) > 0)
            xDBLADD(Q, P, PQdiff, A24p, C24);
        else
            xDBLADD(Q, PQdiff, P, A24p, C24);
    }

    mpz_clear(r);
    mpz_clear(n);
}

void j_invariant(fp2_t j_inv, const fp2_t A, const fp2_t C) {

    // TODO: compare this approach with trying to stay in projective coords
    // until last division

    // Use jinv as register until last line:
    // jinv = A/C: a
    fp2_div_unsafe(j_inv, A, C);

    // t0 = jinv^2: a^2
    fp2_sq_unsafe(lr1, j_inv);

    // jinv = t0 - 3: a^2 - 3
    fp2_sub_uint(j_inv, lr1, 3);

    // t0 = jinv^3: (a^2 - 3)^3
    fp2_sq_unsafe(lr1, j_inv);
    fp2_mul_safe(lr1, j_inv);

    // t1 = jinv - 1: a^2 - 4
    fp2_sub_uint(lr2, j_inv, 1);
    assert(!fp2_is_zero(lr2) && "A was equal 2 or -2");

    // jinv = t0 / t1: (a^2 - 3)^3/(a^2 - 4)
    fp2_div_unsafe(j_inv, lr1, lr2);

    // jinv = 256 * jinv: 256(a^2 - 3)^3/(a^2 - 4)
    // fp2_mul_int(j_inv, j_inv, 256);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
    fp2_add(j_inv, j_inv, j_inv);
}
