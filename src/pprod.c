#include <assert.h>
#include <stdlib.h>

#include "pprod.h"

// Allocate the memory for pprod type and initialize the values
void pprod_init(pprod_t *pp) {
    // Allocate memory for the whole object
    *pp = (pprod_t)malloc(sizeof(struct pprod));

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

void pprod_set_array(pprod_t pp, unsigned int *primes, unsigned int n_primes) {
    if (pp->primes != NULL) {
        free(pp->primes);
        pp->primes = NULL;
    }

    pp->n_primes = n_primes;

    // Allocate the list of numbers
    pp->primes = (unsigned int *)malloc(n_primes * sizeof(unsigned int));

    // Allocate the result of prime primes multiplication
    mpz_init_set_ui(pp->value, 1);

    // Copy all values and calculate the result product
    for (unsigned int i = 0; i < n_primes; i++) {
        assert(primes[i] > 0 && "Given numbers must be larger than 0");

        pp->primes[i] = primes[i];

        // Only odd primes are allowed except the first argument being power of
        // 2 or odd number
        if (i == 0) {
            int bits = 0;
            int p = primes[i];
            while (p) {
                if (p & 1)
                    bits++;
                p /= 2;
            }
            assert((bits == 1 || primes[i] % 2 == 1) &&
                   "First number can only be power of 2 or odd prime");
        }
        if (i > 0)
            assert(primes[i] % 2 && "Only first number can be power of 2, "
                                    "further have to be odd primes");

        mpz_mul_ui(pp->value, pp->value, primes[i]);
    }
}

void pprod_set(pprod_t pp, pprod_t other) {
    pprod_set_array(pp, other->primes, other->n_primes);
}
