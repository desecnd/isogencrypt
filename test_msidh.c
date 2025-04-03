#include <stdlib.h>

#include "testing.h"
#include "proto_msidh.h"

void test_pprod_init() {
    // product: 223'092'870 < 2^32
    unsigned int primes[9] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    pprod_t M;
    pprod_init(&M);

    pprod_set(M, primes, 9);
    CHECK_MSG(!mpz_cmp_ui(M->value, 223092870), "Incorrect product of primes value");

    pprod_clear(&M);
}

void test_random_unit_sampling_large() {
    srand(0xdeafbeef);

    // 100 prime numbers
    unsigned int primes_large[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 
        83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 
        199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 
        347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 
        479, 487, 491, 499, 503, 509, 521, 523, 541
    };

    pprod_t M;
    pprod_init(&M);
    pprod_set(M, primes_large, sizeof(primes_large)/sizeof(unsigned int));

    mpz_t result;
    mpz_init(result);

    // Test if obtained result is really a quadratic root of unity
    for (int i = 0; i < 10; i++) {
        int ret = sample_quadratic_root_of_unity(result, M);
        CHECK_MSG(!ret, "sample_quadratic_root_of_unity returned non-zero value");
        if (ret) break;

        // Check if r^2 == 1 (mod M)
        mpz_mul(result, result, result);
        mpz_mod(result, result, M->value);
        
        CHECK_MSG(!mpz_cmp_ui(result, 1), "Sampling result is not a square root of unity");
    }

    mpz_clear(result);
    pprod_clear(&M);
}

void test_random_unit_sampling_small() {
    // seed random number generators
    // TODO: make sure still using rand in implementation
    srand(0xdeafbeef);

    unsigned int primes_small[5] = {2, 3, 5, 7, 11};
    pprod_t M;
    pprod_init(&M);
    pprod_set(M, primes_small, sizeof(primes_small)/sizeof(unsigned int));

    mpz_t result;
    mpz_init(result);

    // Test if obtained result is really a quadratic root of unity
    for (int i = 0; i < 10; i++) {
        int ret = sample_quadratic_root_of_unity(result, M);
        CHECK_MSG(!ret, "sample_quadratic_root_of_unity returned non-zero value");
        if (ret) break;

        // Check if r^2 == 1 (mod M)
        mpz_mul(result, result, result);
        mpz_mod(result, result, M->value);
        
        CHECK_MSG(!mpz_cmp_ui(result, 1), "Sampling result is not a square root of unity");
    }

    mpz_clear(result);
    pprod_clear(&M);
}

void test_params_generation() {
    mpz_t p;
    mpz_init(p);
    pprod_t A, B;
    pprod_init(&A);
    pprod_init(&B);

    msidh_gen_params(2, p, A, B);
    // 839 is only when 4 is falsely included in Alice primes, otherwise 419 should be used
    CHECK_MSG(!mpz_cmp_ui(p, 839), "Incorrect cofactor and prime number for given msidh params");

    pprod_clear(&A);
    pprod_clear(&B);
    mpz_clear(p);
}

int main() {

    TEST_RUN(test_pprod_init());
    TEST_RUN(test_random_unit_sampling_small());
    TEST_RUN(test_random_unit_sampling_large());
    TEST_RUN(test_params_generation());
    TEST_RUNS_END;

    return 0;
}