#pragma once

#include <stdint.h>

#include "ec_point_xz.h"
#include "fp2.h"

/*
 * @brief Calculate x-coordinate of the double point x(R) = x([2]P).
 * Function is argument-safe for R = P.
 */
void xDBL(point_t R, const point_t P, const fp2_t A24p, const fp2_t C24);

/*
 * @brief Calculate x-coordinate of the point multiplied by the power of 2 x(R)
 * = x([2^e]P). Function is argument-safe for R = P.
 */
void xDBLe(point_t R, const point_t P, const fp2_t A24p, const fp2_t C24,
           const int e);

// Calculate P + Q given P, Q, P - Q
void xADD(point_t PQsum, const point_t P, const point_t Q,
          const point_t PQdiff);

/*
 * @brief Calculate x coordinate of R = [m]P using Montgomery Ladder algorithm.
 *  Function is not argsafe for R0 = P.
 * @ref https://eprint.iacr.org/2017/212.pdf
 */
void xLADDER(point_t R0, const point_t P, const mpz_t m, const fp2_t A24p,
             const fp2_t C24);

/*
 * @brief Calculate x coordinate of R = [m]P using Montgomery Ladder algorithm
 * with m fitting in 63 bits. Function is not argsafe for R0 = P
 * @ref https://eprint.iacr.org/2017/212.pdf
 */
void xLADDER_int(point_t R0, const point_t P, long int m, const fp2_t A24p,
                 const fp2_t C24);

void xLADDER3PT_int(point_t P, point_t Q, point_t PQdiff, long int m,
                    const fp2_t A24p, const fp2_t C24);

void xLADDER3PT(point_t P, point_t Q, point_t PQdiff, const mpz_t m,
                const fp2_t A24p, const fp2_t C24);

/*
 * @brief Calculate j-invariant of the Elliptic Curve in Montgomery Model with
 * coefficient a = (A : C)
 */
void j_invariant(fp2_t j_inv, const fp2_t A, const fp2_t C);

