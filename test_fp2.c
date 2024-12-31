#include <stdio.h>
#include <assert.h>

#include "fp2.h"

int main() {
    printf("[*] Running Fp^2 arithmetic test for char p = 431\n");

    // init characteristic
    fp_t p;
    fp_init(p);
    fp_set_uint(p, 431);
    fp_set_global_char(p);

    fp2_t x, y;
    fp2_init(&x);
    fp2_init(&y);

    // x = 27 + 59i
    fp_set_uint(x->a, 21);
    fp_set_uint(x->b, 59);
    printf("x: %ld, %ldi\n", mpz_get_ui(x->a), mpz_get_ui(x->b));

    // y = 47 + 25i
    fp_set_uint(y->a, 47);
    fp_set_uint(y->b, 25);
    printf("y: %ld, %ldi\n", mpz_get_ui(y->a), mpz_get_ui(y->b));

    fp2_t r; 
    fp2_init(&r);

    // r <- x + y
    // r = (27 + 59i) + (47 + 25i)
    // r = (68 + 84i)
    // r = 68 + 84i (mod 431)
    fp2_add(r, x, y);
    printf("x + y: %ld, %ldi\n", mpz_get_ui(r->a), mpz_get_ui(r->b));
    assert(!mpz_cmp_ui(r->a, 68));
    assert(!mpz_cmp_ui(r->b, 84));

    // r <- x - y
    // r = (27 + 59i) - (47 + 25i)
    // r = (-26 + 34i) 
    // r = 405 + 34i (mod 431)
    fp2_sub(r, x, y);
    printf("x - y: %ld, %ldi\n", mpz_get_ui(r->a), mpz_get_ui(r->b));
    assert(!mpz_cmp_ui(r->a, 405));
    assert(!mpz_cmp_ui(r->b, 34));

    // r <- x * y
    // r = (27 + 59i) * (47 + 25i)
    // r = (-488 + 3298i)
    // r = 374 + 281i (mod 431)
    fp2_mul(r, x, y);
    printf("x * y: %ld, %ldi\n", mpz_get_ui(r->a), mpz_get_ui(r->b));
    assert(!mpz_cmp_ui(r->a, 374));
    assert(!mpz_cmp_ui(r->b, 281));

    // r <- x^2
    // r = (27 + 59i) * (27 + 59i)
    // r = (-3040 + 2478)
    // r = 408 + 323i (mod 431)
    fp2_sq(r, x);
    printf("x ^ 2: %ld, %ldi\n", mpz_get_ui(r->a), mpz_get_ui(r->b));
    assert(!mpz_cmp_ui(r->a, 408));
    assert(!mpz_cmp_ui(r->b, 323));

    fp2_clear(&x);
    fp2_clear(&y);
    fp2_clear(&r);

    printf("[+] All tests passed\n");

    return 0;
}