#pragma once

#include <gmp.h>

#include "ec_tors_basis.h"
#include "pprod.h"

// Number of prime numbers used by a single party
#define TERSIDH_TMIN 2
// TerSIDH t parameter for 128, 192 and 256 bits of security
#define TERSIDH_T128 93
#define TERSIDH_T192 128
#define TERSIDH_T256 162
#define TERSIDH_TMAX 200

// Used by tersidh_state structure
enum {
    TERSIDH_STATUS_UNINITIALIZED = 0,
    TERSIDH_STATUS_INITIALIZED,
    TERSIDH_STATUS_PREPARED,
    TERSIDH_STATUS_EXCHANGED
};

struct tersidh_state {
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
    // if secret is set to nonzero it will be 
    mpz_t secret;
    point_t KP, KQ;
    pprod_t KP_deg, KQ_deg;

    // Shared Key -> j_invariant
    fp2_t j_inv;

    // Current state of the protocol
    int status;
};

struct tersidh_data {
    int t, f;
    fp2_t a, xP, xQ, xR;
};

struct tersidh_const_data {
    int t, f;
    const char *a_str, *xP_str, *xQ_str, *xR_str;
};

void tersidh_data_copy_const(struct tersidh_data *td, const struct tersidh_const_data *tcd);

void tersidh_data_init(struct tersidh_data *td);

void tersidh_data_clear(struct tersidh_data *td);

void tersidh_state_init(struct tersidh_state *tersidh);

void tersidh_state_clear(struct tersidh_state *tersidh);

void tersidh_state_reset(struct tersidh_state *tersidh);

void tersidh_state_prepare(struct tersidh_state *tersidh,
                         const struct tersidh_data *params, int is_bob);

void tersidh_key_exchange(struct tersidh_state *tersidh,
                        const struct tersidh_data *pk_other);

void tersidh_get_pubkey(const struct tersidh_state *tersidh,
                      struct tersidh_data *pk_self);

void tersidh_generate_kernel_points(struct tersidh_state* tersidh, int skip_secret);

/*
 * @brief Given security parameter t, and cofactor generate public params used
 * in TERSIDH: p, A, B, where p = fAB - 1 is prime. Return cofactor f or -1 if
 * something gone wrong. This function searches for the cofactor in range. Prime
 * p is always congruent to 3 mod 4 due to the construction of the prime (4 | A
 * => 4 | ABf)
 */
int tersidh_gen_pub_params(mpz_t p, pprod_t A, pprod_t B, int t);

/*
 * @brief Given security parameter t and cofactor f, calculate public params
 */
int tersidh_calc_pub_params(mpz_t p, pprod_t A, pprod_t B, int t, int f);

/*
 * @brief Generate TERSIDH public key from Alice perspective
 */

void _tersidh_gen_pubkey_alice(fp2_t A24p_alice, fp2_t C24_alice, point_t KP, point_t KQ,
                             struct tors_basis *PQ_bob, const pprod_t A_deg, const fp2_t A24p_base,
                             const fp2_t C24_base);
// void tersidh_gen_pubkey(fp2_t A24p_alice, fp2_t C24_alice, struct tors_basis*
// PQ_alice, struct tors_basis* PQ_bob, const fp2_t A24p_base, const fp2_t
// C24_base, const mpz_t secret, const mpz_t mask);

/*
 * @brief Calculate shared secret (j_inv) and destination curve coefficient (a +
 * 2)/4 = (A24p : C24) in xDBL form
 */
void _tersidh_key_exchange_alice(fp2_t j_inv, fp2_t A24p_final, fp2_t C24_final,
                               const fp2_t A24p_bob, const fp2_t C24_bob,
                               struct tors_basis *BPQA, const pprod_t A, const mpz_t A_sec);
