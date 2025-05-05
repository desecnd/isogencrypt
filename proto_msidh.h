#ifndef PROTO_MSIDH_H
#define PROTO_MSIDH_H

#include <gmp.h>

#include "pprod.h"
#include "ec_mont.h"

struct tors_basis {
    point_t P, Q, PQd;
    mpz_t n;
};

int sample_quadratic_root_of_unity(mpz_t result, pprod_t modulus);

/*
 * @brief Given security parameter t, generate public params used in MSIDH: p, A, B, where p = fAB - 1 is prime.
 *  Return cofactor f or -1 if something gone wrong.
 */
int msidh_gen_pub_params(mpz_t p, pprod_t A, pprod_t B, unsigned int t);

/* 
 * @brief Generate MSIDH public key from Alice perspective
*/
void msidh_genkey(fp2_t A24p_alice, fp2_t C24_alice, point_t PB, point_t QB, point_t PQBdiff, const mpz_t p, const fp2_t A24p_base, const fp2_t C24_base, point_t PA, point_t QA, point_t PQAdiff, const pprod_t deg_A, const pprod_t deg_B, const mpz_t secret, const mpz_t mask);

/*
 * @brief Calculate subgroup basis of the torsion basis (Pn, Qn) = [N/n](P, Q) of order n, where N is the order of (P, Q)
 */
void tors_basis_get_subgroup(struct tors_basis *PQsub, mpz_t n, const struct tors_basis *PQ, const fp2_t A24p, const fp2_t C24);

#endif