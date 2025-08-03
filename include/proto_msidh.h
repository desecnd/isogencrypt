#pragma once

#include <gmp.h>

#include "ec_tors_basis.h"
#include "pprod.h"

// Used by msidh_state structure
enum {
    MSIDH_STATUS_UNINITIALIZED = 0,
    MSIDH_STATUS_INITIALIZED,
    MSIDH_STATUS_PREPARED,
    MSIDH_STATUS_EXCHANGED
};

struct msidh_state {
    gmp_randstate_t randstate;

    int t, f;
    int is_bob;

    mpz_t p;
    pprod_t A, B;

    // Torsion basis for Alice and for Bob
    struct tors_basis PQ_self, PQ_pubkey;

    // Starting Elliptic Curve coefficient a = A/C in xDBL form
    fp2_t A24p_start, C24_start;

    // Public Key Elliptic Curve
    fp2_t A24p_pubkey, C24_pubkey;

    // Private Key
    mpz_t secret;

    // Shared Key -> j_invariant
    fp2_t j_inv;

    // Current state of the protocol
    int status;
};

struct msidh_data {
    int t, f;
    fp2_t a, xP, xQ, xR;
};

struct msidh_const_data {
    int t, f;
    const char *a_str, *xP_str, *xQ_str, *xR_str;
};

// TODO: From MSIDH paper: 128bit prime: p = 2^2 * l1 ... l571 * 10 - 1
// So at least 571//2 primes per side is required
#define MSIDH_TMIN 4
#define MSIDH_TMAX 600
#define MSIDH_TMAX_HALF MSIDH_TMAX / 2

void msidh_data_init(struct msidh_data *md);

void msidh_data_clear(struct msidh_data *md);

void msidh_state_init(struct msidh_state *msidh);

void msidh_state_clear(struct msidh_state *msidh);

void msidh_state_reset(struct msidh_state *msidh);

void msidh_state_prepare(struct msidh_state *msidh,
                         const struct msidh_data *params, int is_bob);

void msidh_key_exchange(struct msidh_state *msidh,
                        const struct msidh_data *pk_other);

void msidh_get_pubkey(const struct msidh_state *msidh,
                      struct msidh_data *pk_self);

// TODO: add gmp_randstate instead of random
int sample_quadratic_root_of_unity(mpz_t result, pprod_t modulus);

/*
 * @brief Given security parameter t, and cofactor generate public params used in MSIDH: p,
 * A, B, where p = fAB - 1 is prime. Return cofactor f or -1 if something gone
 * wrong. This function searches for the cofactor in range. Prime p is always congruent to 3 mod 4 due to the construction of the
 * prime (4 | A => 4 | ABf)
 */
int msidh_gen_pub_params(mpz_t p, pprod_t A, pprod_t B, int t);

/* 
 * @brief Given security parameter t and cofactor f, calculate public params
 */
int msidh_calc_pub_params(mpz_t p, pprod_t A, pprod_t B, int t, int f);

/*
 * @brief Generate MSIDH public key from Alice perspective
 */
void _msidh_gen_pubkey_alice(fp2_t A24p_alice, fp2_t C24_alice,
                             struct tors_basis *PQ_alice,
                             struct tors_basis *PQ_bob, const fp2_t A24p_base,
                             const fp2_t C24_base, const mpz_t secret,
                             const mpz_t mask);
// void msidh_gen_pubkey(fp2_t A24p_alice, fp2_t C24_alice, struct tors_basis*
// PQ_alice, struct tors_basis* PQ_bob, const fp2_t A24p_base, const fp2_t
// C24_base, const mpz_t secret, const mpz_t mask);

/*
 * @brief Calculate shared secret (j_inv) and destination curve coefficient (a +
 * 2)/4 = (A24p : C24) in xDBL form
 */
void _msidh_key_exchange_alice(fp2_t j_inv, fp2_t A24p_final, fp2_t C24_final,
                               const fp2_t A24p_bob, const fp2_t C24_bob,
                               struct tors_basis *BPQA, const mpz_t A_sec);
