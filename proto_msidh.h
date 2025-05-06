#ifndef PROTO_MSIDH_H
#define PROTO_MSIDH_H

#include <gmp.h>

#include "pprod.h"
#include "ec_mont.h"

struct tors_basis {
    point_t P, Q, PQd;
    pprod_t n;
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
void msidh_genkey(fp2_t A24p_alice, fp2_t C24_alice, point_t PB, point_t QB, point_t PQBdiff, const fp2_t A24p_base, const fp2_t C24_base, point_t PA, point_t QA, point_t PQAdiff, const pprod_t deg_A, const pprod_t deg_B, const mpz_t secret, const mpz_t mask);

/*
 * @brief Calculate shared secret (j_inv) and destination curve coefficient (a + 2)/4 = (A24p : C24) in xDBL form
 */
void msidh_key_exchange(fp2_t j_inv, fp2_t A24p_final, fp2_t C24_final, const fp2_t A24p_bob, const fp2_t C24_bob, const struct tors_basis* BPQA, const mpz_t A_sec); 

/*
 * @brief Calculate subgroup basis of the torsion basis (Pn, Qn) = [N/n](P, Q) of order n, where N is the order of (P, Q)
 */
void tors_basis_get_subgroup(struct tors_basis *PQsub, pprod_t n, const struct tors_basis *PQ, const fp2_t A24p, const fp2_t C24);

#endif