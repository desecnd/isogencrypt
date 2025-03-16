#ifndef EC_MONT_H
#define EC_MONT_H

#include "fp2.h"

// Point in projective coordinates X-Z
struct point_xz { 
    fp2_t X; 
    fp2_t Z; 
};

typedef struct point_xz* point_t;

void point_set(point_t R, const point_t P);

void point_init(point_t *P);

void point_clear(point_t *P);

// Set: x(P) = x, and z(P) = 1
void point_set_str_x(point_t P, const char *x); 

// Normalize: P = (X : Z) -> (X' : 1) with X' = X/Z
void point_normalize_coords(point_t P);

// Set A24p := A + 2C and C24 := 4C
static inline void calc_curve_proj_coeffs(fp2_t A24_plus, fp2_t C24, const fp2_t A, const fp2_t C) {
    fp2_set(A24_plus, A);               // A24p = A
    fp2_add(C24, C, C);                 // C24 = 2C
    fp2_add(A24_plus, A24_plus, C24);   // A24p = A + 2C
    fp2_add(C24, C24, C24);             // C24 = 4C
}

// Calculate the double of the point
void xDBL(point_t R, const point_t P, const fp2_t A24_plus, const fp2_t C24);

// Calculate R: multiple [2^e]P of point P
void xDBLe(point_t R, const point_t P, const fp2_t A24_plus, const fp2_t C24, const int e);

void criss_cross(fp2_t lsum, fp2_t rdiff, const fp2_t x, const fp2_t y, const fp2_t z, const fp2_t w);

// Calculate P + Q given P, Q, P - Q
void xADD(point_t PQsum, const point_t P, const point_t Q, const point_t PQdiff);

void xLADDER3PT_int(point_t P, point_t Q, point_t PQdiff, long int m, const fp2_t A24p, const fp2_t C24);
void xLADDER3PT(point_t P, point_t Q, point_t PQdiff, mpz_t m, const fp2_t A24p, const fp2_t C24);



#endif