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

/* 
 * @brief Transform (A : C) = (a : 1) into (A' : C') = ((a+ 2)/4 : 1) = (a + 2 : 4) suitable for xDBL
 */
void A24p_from_A(fp2_t A24_plus, fp2_t C24, const fp2_t A, const fp2_t C);

/* 
 * @brief Transform (A' : C') = (a + 2 : 4) into (A : C) = (a : 1) 
 */
void A_from_A24p(fp2_t A, fp2_t C, const fp2_t A24_plus, const fp2_t C24);

// Calculate the double of the point
void xDBL(point_t R, const point_t P, const fp2_t A24_plus, const fp2_t C24);

// Calculate R: multiple [2^e]P of point P
void xDBLe(point_t R, const point_t P, const fp2_t A24_plus, const fp2_t C24, const int e);

void criss_cross(fp2_t lsum, fp2_t rdiff, const fp2_t x, const fp2_t y, const fp2_t z, const fp2_t w);

// Calculate P + Q given P, Q, P - Q
void xADD(point_t PQsum, const point_t P, const point_t Q, const point_t PQdiff);

void xLADDER3PT_int(point_t P, point_t Q, point_t PQdiff, long int m, const fp2_t A24p, const fp2_t C24);

void xLADDER3PT(point_t P, point_t Q, point_t PQdiff, mpz_t m, const fp2_t A24p, const fp2_t C24);

/*
 * @brief Return size of the kernel points array based on odd isogeny degree
 */
#define KPS_DEG2SIZE(__deg) (__deg - 1)/2

/*
 * @brief Generate points that lay in the kernel of the isogeny of degree d using generator P
 * @param[out]  kpoints     list of `P` multiples up to `deg`: `[[1]P, [2]P, [3]P, ... [d - 1]P]`
 * @param[in]   n           size of the kpoints array (degree of the isogeny - 1)
 * @see https://eprint.iacr.org/2017/504
*/
void KPS(point_t * kpoints, size_t n, const point_t P, const fp2_t A24p, const fp2_t C24);

/*
 * @brief Transform (X : Z) kernel points into (X + Z : X - Z) pairs for efficient computation
*/
void prepare_kernel_points(point_t *kpts, size_t n);

/*
 * @brief Calculate x-coordinate of the odd-degree isogeny point image Q = φ(P) given ker(phi) = prep_points
 */
void xISOG_point(point_t Q, const point_t *prep_kpts, size_t n, const point_t P);


/*
 * @brief Calculate a coefficient of the odd-degree isogeny curve codomain E' = φ(E) given kernel points list
 */
void aISOG_curve_KPS(fp2_t A_, fp2_t C_, const fp2_t A24p, const fp2_t C24, const point_t * kpts, size_t n);

/*
 * @brief Calculate a coefficient of the odd-degree isogeny curve codomain E' = φ(E) given kernel point
 * @details
 *  This function allocates space and generates KPS list, later calling `aISOG_curve_KPS` to obtain the result.
 *  If KPS is calulated beforehand use `aISOG_curve_KPS` istead to avoid unnecessary computations.
 */
void aISOG_curve(fp2_t A_, fp2_t C_, const fp2_t A24p, const fp2_t C24, const point_t K, int degree);

#endif