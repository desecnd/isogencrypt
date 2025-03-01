#ifndef PROTO_MSIDH_H
#define PROTO_MSIDH_H

#include <gmp.h>

// Structure representing smooth composite number being product of small primes.
// In factorization each exponent is equal to 1, so each factor
// should be present only once
//
// in MSDIH max prime factor for 256-bit variant 
// would be around 1600
struct pprod {
    mpz_t value;
    unsigned int *primes;
    unsigned int n_primes;
};
typedef struct pprod *pprod_t;

void pprod_init(pprod_t *pp, unsigned int * primes, unsigned int n_primes);
void pprod_clear(pprod_t *pp);

int sample_quadratic_root_of_unity(mpz_t result, pprod_t modulus);

#endif