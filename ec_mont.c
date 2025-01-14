#include <stdlib.h>
#include <string.h>

#include "ec_mont.h"

// copy coordinates of point: R <- P
void point_set(point_t R, const point_t P) {
    // R.X = P.X, R.Y = P.Y
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

// Initialize the point with zP = 1 and xP = xa + xb * i
void point_set_str(point_t P, const char *xa, const char *xb) {
    fp2_set_str(P->X, xa, xb);
    fp2_set_uint(P->Z, 1);
}


// Calculate the double of the point
void xDBL(point_t R, const point_t P, const fp2_t A24_plus, const fp2_t C24) {
    // By X, Z we mean x(P) and z(P)
    fp2_t t0, t1;
    fp2_init(&t0);
    fp2_init(&t1);

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
    // Function is argument safe for calling PQSum = Q or P
    // Given P, Q and P-Q output P+Q
    // Cost: 4M + 2S + 6a

    fp2_t t0, t1;
    fp2_init(&t0); fp2_init(&t1);

    fp2_add(t0, P->X, P->Z); // t0 = xP + zP
    fp2_sub(t1, P->X, P->Z); // t1 = xP - zP
    fp2_sub(PQsum->X, Q->X, Q->Z); // xP = xQ - zQ
    fp2_add(PQsum->Z, Q->X, Q->Z); // zP = xQ + zQ

    fp2_mul_safe(t0, PQsum->X); // t0 *= xP
    fp2_mul_safe(t1, PQsum->Z); // t1 *= zP

    fp2_sub(PQsum->Z, t0, t1); // zP = t0 - t1
    fp2_add(PQsum->X, t0, t1); // xP = t0 + t1

    fp2_sq_unsafe(PQsum->Z, P->Z);
    fp2_sq_unsafe(PQsum->X, P->X);// XQP = XP * XP

    fp2_mul_safe(PQsum->Z, PQdiff->X); // ZQP = xPQ * ZP
    fp2_mul_safe(PQsum->X, PQdiff->Z); // XQP = XQP * zPQ

    fp2_clear(&t0); fp2_clear(&t1);
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

// TODO: substitute with unlimited integer later
// calculate P = P + [m]Q
void LADDER3PT_uint(point_t P, point_t Q, point_t PQdiff, long int m, const fp2_t A24p, const fp2_t C24) {

    // For X-only Montgomery arithmetic [-m]P = [m]P
    // as -P only influences y-coordinate
    m = (m < 0) ? -m : m;

    while (m > 0) {
        if (m & 1) xDBLADD(Q, P, PQdiff, A24p, C24);
        else xDBLADD(Q, PQdiff, P, A24p, C24);
    }

}