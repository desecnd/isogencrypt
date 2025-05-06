#ifndef PROTO_MSIDH_H
#define PROTO_MSIDH_H

#include <gmp.h>

#include "pprod.h"
#include "ec_mont.h"

struct tors_basis {
    point_t P, Q, PQd;
    pprod_t n;
};

// Used by msidh_state structure
enum {
    MSIDH_STATUS_UNINITIALIZED = 0,
    MSIDH_STATUS_INITIALIZED,
    MSIDH_STATUS_PREPARED, 
    MSIDH_STATUS_EXCHANGED 
};

struct msidh_state {
    gmp_randstate_t randstate;

    struct tors_basis PQ_me, PQ_other;
    // Public Key
    fp2_t pk_A24p, pk_C24;
    // Private Key
    mpz_t secret;
    // Shared Key -> j_invariant
    fp2_t sk_jinv;

    int status; 
};

// TODO: Make sure that the global arithmetic is set to FP
void msidh_state_init(struct msidh_state* msidh);

void msidh_state_clear(struct msidh_state* msidh);

void msidh_prepare(struct msidh_state *msidh, const fp2_t A24p, const fp2_t C24, const struct tors_basis *PQ, pprod_t A_deg, pprod_t B_deg, int is_bob);

// TODO: add gmp_randstate instead of random
int sample_quadratic_root_of_unity(mpz_t result, pprod_t modulus);

/*
 * @brief Given security parameter t, generate public params used in MSIDH: p, A, B, where p = fAB - 1 is prime.
 *  Return cofactor f or -1 if something gone wrong.
 */
int msidh_gen_pub_params(mpz_t p, pprod_t A, pprod_t B, unsigned int t);

/* 
 * @brief Generate MSIDH public key from Alice perspective
*/
void _msidh_gen_pubkey_alice(fp2_t A24p_alice, fp2_t C24_alice, struct tors_basis* PQ_alice, struct tors_basis* PQ_bob, const fp2_t A24p_base, const fp2_t C24_base, const mpz_t secret, const mpz_t mask);
// void msidh_gen_pubkey(fp2_t A24p_alice, fp2_t C24_alice, struct tors_basis* PQ_alice, struct tors_basis* PQ_bob, const fp2_t A24p_base, const fp2_t C24_base, const mpz_t secret, const mpz_t mask);

/*
 * @brief Calculate shared secret (j_inv) and destination curve coefficient (a + 2)/4 = (A24p : C24) in xDBL form
 */
void _msidh_key_exchange_alice(fp2_t j_inv, fp2_t A24p_final, fp2_t C24_final, const fp2_t A24p_bob, const fp2_t C24_bob, struct tors_basis * BPQA, const mpz_t A_sec);

void msidh_key_exchange(struct msidh_state *msidh, const fp2_t A24p, const fp2_t C24, const struct tors_basis* BPQA);

/*
 * @brief Calculate subgroup basis of the torsion basis (Pn, Qn) = [N/n](P, Q) of order n, where N is the order of (P, Q)
 */
void tors_basis_get_subgroup(struct tors_basis *PQsub, pprod_t n, const struct tors_basis *PQ, const fp2_t A24p, const fp2_t C24);

#endif