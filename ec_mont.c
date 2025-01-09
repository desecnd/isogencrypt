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

// Calculate the double of the point
void xDBL(point_t R, const point_t P, const fp2_t A24_plus, const fp2_t C24) {
    // By X, Z we mean x(P) and z(P)
    fp2_t t0, t1;
    fp2_init(&t0);
    fp2_init(&t1);

    fp2_sub(t0, P->X, P->Z);    // t0 = X - Z
    fp2_add(t1, P->X, P->Z);    // t1 = X + Z

    // Warning! Must use "safe" function for square-ing
    fp2_sq_safe(t0, t0);             // t0 = (X - Z)^2 
    fp2_sq_safe(t1, t1);             // t1 = (X + Z)^2

    fp2_mul_unsafe(R->Z, t0, C24);     // Z' = (X - Z)^2 * C24
    fp2_mul_unsafe(R->X, R->Z, t1);    // X' = (X + Z)^2 * (X - Z)^2 * C24

    fp2_sub(t1, t1, t0);        // t1 = (X + Z)^2 - (X - Z)^2 
    fp2_mul_unsafe(t0, A24_plus, t1);  // t0 = A24p * [(X + Z)^2 - (X - Z)^2]

    fp2_add(R->Z, R->Z, t0);    // Z' = A24p * [(X + Z)^2 - (X - Z)^2] + (X - Z)^2 * C24
    fp2_mul_safe(R->Z, R->Z, t1);    // Z' = { A24p * [(X + Z)^2 - (X - Z)^2] + C24 (X - Z)^2] } [(X + Z)^2 - (X - Z)^2]

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


// Calculate the codomain (A24+, C24) of the 2-degree isogeny from given kernel K of order 2
void isog2_codomain(const point_t K, fp2_t A24_plus, fp2_t C24) {
    // A24+ = x(K)^2
    fp2_sq_unsafe(A24_plus, K->X);
    // C24 = z(K)^2
    fp2_sq_unsafe(C24, K->Z);
    // A24+ = z(K)^2 - x(K)^2
    fp2_sub(A24_plus, C24, A24_plus);
}

// Calculate x(P + kQ) given x(P + Q) and z(P + Q)
// void xDBLADD(point_xz P, point_xz Q, const fp2_t XPQ, const fp2_t ZPQ, const fp2_t A24) {
