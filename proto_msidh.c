#include <stdlib.h>
#include <gmp.h>
#include <stdio.h>

#include "fp2.h"
#include "ec_mont.h"
#include "proto_msidh.h"

// Allocate the memory for pprod type and initialize the values
void pprod_init(pprod_t *pp, unsigned int * primes, unsigned int n_primes) {
    // Allocate memory for the whole object
    *pp = (pprod_t) malloc(sizeof(struct pprod));

    // Allocate the list of numbers
    (*pp)->primes = (unsigned int*) malloc(n_primes * sizeof(unsigned int));
    // Allocate the result of prime primes multiplication
    mpz_init_set_ui((*pp)->value, 1);

    // Copy all values and calculate the result product
    for (unsigned int i = 0; i < n_primes; i++) {
        (*pp)->primes[i] = primes[i];
        mpz_mul_ui((*pp)->value, (*pp)->value, primes[i]);
    }
    (*pp)->n_primes = n_primes;
}

// Deallocate the memory for pprod type
void pprod_clear(pprod_t *pp) {
    mpz_clear((*pp)->value);
    free((*pp)->primes);
    free(*pp);
    pp = NULL;
}


// Sample an element `x` from ``Z/mZ`` where ``x^2 = 1 (mod m)``.
// Return != 0 if something goes wrong
int sample_quadratic_root_of_unity(mpz_t result, pprod_t modulus) {
    // TODO: CRT: check if it will overflow?
    // Use fp_t instead of mpz_t if so

    int ret = 0;
    mpz_set_ui(result, 0);

    mpz_t m, inv;
    mpz_init(m);
    mpz_init(inv);

    for (unsigned i = 0; i < modulus->n_primes; i++) {
        unsigned int p = modulus->primes[i];

        // Toss a coin and choose 1 or -1 (== p - 1)
        // TODO: change rand to something more serious in the future
        unsigned int root = 1;
        if (rand() % 2 == 0) {
            root = p - 1;
        }

        // Solve CRT Congruence: root^2 == 1 (mod p)

        // m = M / p 
        mpz_div_ui(m, modulus->value, p);

        // set temporary: inv = p
        // inv = (M/p)^-1 (mod p)
        // if returns 0: inverse does not exist
        mpz_set_ui(inv, p);
        int inv_ret = mpz_invert(inv, m, inv);
        if (inv_ret == 0) {
            ret = -1;
            break;
        }

        // m = root * (M / p)
        mpz_mul_ui(m, m, root);

        // result += root * (M/p) * (M/p)^-1
        mpz_addmul(result, m, inv);
        mpz_mod(result, result, modulus->value);
    }
    mpz_clear(m);
    mpz_clear(inv);

    return ret;
}



