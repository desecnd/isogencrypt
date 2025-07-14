#pragma once

#include <gmp.h>

/*
 * @class Structure representing smooth composite number being product of small
 * primes.
 * @brief In factorization each exponent is equal to 1, so each factor
 * should be present only once
 *
 * in MSDIH max prime factor for 256-bit variant
 * would be around 1600
 */
struct pprod {
    mpz_t value;
    unsigned int *primes;
    unsigned int n_primes;
};
typedef struct pprod *pprod_t;

void pprod_init(pprod_t *pp);
void pprod_clear(pprod_t *pp);

/*
 * @brief Set value of pprod_t number by taking a product of primes. Only the
 * first argument of the `primes` table can be (optionally) a power of 2.
 */
void pprod_set_array(pprod_t pp, unsigned int *primes, unsigned int n_primes);

/*
 * @brief Set value of pprod_t number by calling pprod_set_array on data of the
 * other pprod_t
 */
void pprod_set(pprod_t pp, pprod_t other);
