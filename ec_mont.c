#include <stdlib.h>
#include <string.h>

#include "ec_mont.h"

// copy coordinates of point: R <- P
void point_set(point_t R, const point_t P) {
    // R.X = P.X, R.Y = P.Y

    // TODO: Not a deepcopy -> only copies the pointers
    memcpy(R, P, sizeof(struct point_xz));
}

void point_init(point_t *P) {
    *P = (point_t) malloc(sizeof(struct point_xz));
    // Recursively allocate memory for the member pointers
    fp2_init(&(*P)->X);
    fp2_init(&(*P)->Z);
}

void point_clear(point_t *P) {
    if (P) {
        fp2_clear(&(*P)->X);
        fp2_clear(&(*P)->Z);
        free(*P);
        P = NULL;
    }
}

// Set: x(P) = x, and z(P) = 1
void point_set_str_x(point_t P, const char *x) {
    fp2_set_str(P->X, x);
    fp2_set_uint(P->Z, 1);
}

// Normalize: P = (X : Z) -> (X' : 1) with X' = X/Z
void point_normalize_coords(point_t P) {
    // Registers: 1

    fp2_t t; fp2_init(&t);

    fp2_div_unsafe(t, P->X, P->Z);
    fp2_set(P->X, t);
    fp2_set_uint(P->Z, 1);

    fp2_clear(&t);
}

// Calculate the double of the point R = [2]P
void xDBL(point_t R, const point_t P, const fp2_t A24_plus, const fp2_t C24) {
    fp2_t t0, t1;
    fp2_init(&t0);
    fp2_init(&t1);

    // P = (X : Z)
    fp2_sub(t0, P->X, P->Z);    // t0 = X - Z
    fp2_add(t1, P->X, P->Z);    // t1 = X + Z

    // Warning! Must use "safe" function for square-ing
    fp2_sq_safe(t0);             // t0 = (X - Z)^2 
    fp2_sq_safe(t1);             // t1 = (X + Z)^2

    fp2_mul_unsafe(R->Z, t0, C24);     // Z' = (X - Z)^2 * C24
    fp2_mul_unsafe(R->X, R->Z, t1);    // X' = (X + Z)^2 * (X - Z)^2 * C24

    fp2_sub(t1, t1, t0);        // t1 = (X + Z)^2 - (X - Z)^2 
    fp2_mul_unsafe(t0, A24_plus, t1);  // t0 = A24p * [(X + Z)^2 - (X - Z)^2]

    fp2_add(R->Z, R->Z, t0);    // Z' = A24p * [(X + Z)^2 - (X - Z)^2] + (X - Z)^2 * C24
    fp2_mul_safe(R->Z, t1);    // Z' = { A24p * [(X + Z)^2 - (X - Z)^2] + C24 (X - Z)^2] } [(X + Z)^2 - (X - Z)^2]

    fp2_clear(&t0);
    fp2_clear(&t1);
}

// Calculate R: multiple [2^e]P of point P
void xDBLe(point_t R, const point_t P, const fp2_t A24_plus, const fp2_t C24, const int e) {

    // Set R <- P
    point_set(R, P);

    // Repeat the step of doubling multiple times
    for (int i = 0; i < e; i++) {
        xDBL(R, R, A24_plus, C24);
    }
}

void xADD(point_t PQsum, const point_t P, const point_t Q, const point_t PQdiff) {
    // Function is argument safe for calling PQSum = Q or P or PQdiff
    // Given P, Q and PQdiff = P-Q output P+Q

    // xADD formula in projective coordinates: P + Q = (X' : Z')
    // X' = ZR * (XP * XQ - ZP * ZQ)^2 
    // Z' = XR * (XP * ZQ - XQ * ZP)^2
    // 
    // Optimised version:
    // X' = ZR * [(XP - ZP)(XQ + ZQ) + (XP + ZP)(XQ - ZQ)]^2 
    // Z' = XR * [(XP - ZP)(XQ + ZQ) - (XP + ZP)(XQ - ZQ)]^2 

    fp2_t t0, t1, t2, t3;
    fp2_init(&t0); fp2_init(&t1); fp2_init(&t2); fp2_init(&t3); 

    // Defining additional variables:
    // a: XP + ZP
    // b: XP - ZP
    // c: XQ + ZQ
    // d: XQ - ZQ
    //
    // X' = ZR * (b * c + a * d)^2 
    // Z' = XR * (b * c - a * d)^2 
    fp2_add(t0, P->X, P->Z);        // t0 = a: XP + ZP
    fp2_sub(t1, P->X, P->Z);        // t1 = b: XP - ZP
    fp2_add(t2, Q->X, Q->Z);        // t2 = c: xQ + zQ
    fp2_sub(t3, Q->X, Q->Z);        // t3 = d: xQ - zQ

    fp2_mul_safe(t0, t3);           // t0 = t0 * t3: a * d
    fp2_mul_safe(t1, t2);           // t1 = t1 * t2: b * c

    fp2_add(t2, t0, t1);            // t2 = t0 + t1: ad + bc
    fp2_sub(t3, t0, t1);            // t3 = t0 - t1: ad - bc

    fp2_sq_unsafe(t0, t2);          // t0 = t2^2: (ad + bc)^2
    fp2_sq_unsafe(t1, t3);          // t1 = t3^2: (ad - bc)^2

    // We must use t2 to not override the used later PQdiff->X in scenario when PQsum = PQdiff 
    // t2 = t0 * Z(P-Q): (ad + bc)^2 * Z(P-Q)
    fp2_mul_unsafe(t2, t0, PQdiff->Z);          
    // Z' = t1 * X(P-Q): (ad - bc)^2 * Z(P-Q)
    fp2_mul_unsafe(PQsum->Z, t1, PQdiff->X);   
    // X' = t2
    fp2_set(PQsum->X, t2);

    fp2_clear(&t0); fp2_clear(&t1); fp2_clear(&t2); fp2_clear(&t3);
}

// TODO: check argument-safeness
// TODO: think about expressing function 
// as set of arithmetic instructions 
// instead of combinatin of 2 function calls
void xDBLADD(point_t P, point_t Q, const point_t PQdiff, const fp2_t A24p, const fp2_t C24) {
    // Out: Q = P + C, Q = 2P
    xADD(Q, P, Q, PQdiff);
    xDBL(P, P, A24p, C24);
}

void criss_cross(fp2_t lsum, fp2_t rdiff, const fp2_t x, const fp2_t y, const fp2_t z, const fp2_t w) {
    // Calculate (xw + yz, xw - yz) given (x, y, z, w)
    // Argument-safe: Yes
    // Registers: 2
    // Cost: 2M + 2a

    fp2_t t0, t1;
    fp2_init(&t0); fp2_init(&t1);

    // t0 = x * w
    fp2_mul_unsafe(t0, x, w);
    // t1 = y * z
    fp2_mul_unsafe(t1, y, z);
    
    // rdiff = t0 - t1: xw - yz
    fp2_sub(rdiff, t0, t1);
    // lsum = t0 + t1: xw + yz
    fp2_add(lsum, t0, t1);
    
    fp2_clear(&t0); fp2_clear(&t1);
}

// calculate P = P + [m]Q
void xLADDER3PT_int(point_t P, point_t Q, point_t PQdiff, long int m, const fp2_t A24p, const fp2_t C24) {
    assert(m > 0 && "Given scalar m must be nonnegative");

    while (m > 0) {
        if (m & 1) 
            xDBLADD(Q, P, PQdiff, A24p, C24);
        else       
            xDBLADD(Q, PQdiff, P, A24p, C24);
        m /= 2;
    }
}

// calculate P = P + [m]Q
void xLADDER3PT(point_t P, point_t Q, point_t PQdiff, mpz_t m, const fp2_t A24p, const fp2_t C24) {
    assert(mpz_sgn(m) > 0 && "Given scalar m must be nonnegative");

    mpz_t r;
    mpz_init(r);

    while (mpz_sgn(m) > 0) {
        // m = m//2; r = m % 2
        mpz_fdiv_qr_ui(m, r, m, 2);

        if (mpz_sgn(r) > 0)
            xDBLADD(Q, P, PQdiff, A24p, C24);
        else 
            xDBLADD(Q, PQdiff, P, A24p, C24);
    }

    mpz_clear(r);
}


