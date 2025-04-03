#include <stdlib.h>
#include "pprod.h"

// Allocate the memory for pprod type and initialize the values
void pprod_init(pprod_t *pp) {
    // Allocate memory for the whole object
    *pp = (pprod_t) malloc(sizeof(struct pprod));

    (*pp)->n_primes = 0;
    (*pp)->primes = NULL;
    mpz_init((*pp)->value);
}


// Deallocate the memory for pprod type
void pprod_clear(pprod_t *pp) {
    mpz_clear((*pp)->value);
    if ((*pp)->primes != NULL) {
        free((*pp)->primes);
    }
    free(*pp);
    pp = NULL;
}

// Set pp = prod(primes)
void pprod_set(pprod_t pp, unsigned int *primes, unsigned int n_primes) {
    if (pp->primes != NULL) {
        free(pp->primes);
        pp->primes = NULL;
    }

    pp->n_primes = n_primes;

    // Allocate the list of numbers
    pp->primes = (unsigned int*) malloc(n_primes * sizeof(unsigned int));

    // Allocate the result of prime primes multiplication
    mpz_init_set_ui(pp->value, 1);

    // Copy all values and calculate the result product
    for (unsigned int i = 0; i < n_primes; i++) {
        pp->primes[i] = primes[i];
        mpz_mul_ui(pp->value, pp->value, primes[i]);
    }
}


