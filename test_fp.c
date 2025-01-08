#include <stdio.h>
#include <assert.h>

#include "fp2.h"

int test_simple() {
    fp_char_setup_uint(431);

    fp_t r;
    fp_init(r);

    fp_set_uint(r, 16); // r == 16
    assert(!mpz_cmp_ui(r, 16));

    fp_sqrt(r, r); // r == 4
    // TODO: incosistent api, change later
    assert(!mpz_cmp_ui(r, 4));

    fp_add_uint(r, r, 10); // r == 14
    assert(!mpz_cmp_ui(r, 14));

    fp_sub_uint(r, r, 15); // r == -1 == 430
    assert(!mpz_cmp_ui(r, 430));

    fp_inv(r, r); // r == r ^-1 == 430
    assert(!mpz_cmp_ui(r, 430));

    fp_clear(r);

    assert(!fp_char_clear());

    return 0;
}

void test_modulo() {
    assert(!fp_char_setup_uint(431));

    fp_t a, b;
    fp_init(a);
    fp_init(b);
    fp_set_uint(a, 512);
    fp_set_uint(b, 791);

    {
        fp_t r;
        fp_init(r);
        // r = a + b (mod 431)
        fp_add(r, a, b);
        assert(!mpz_cmp_ui(r, 10));

        fp_sub(r, a, b);
        assert(!mpz_cmp_ui(r, 152));

        fp_sub(r, b, a);
        assert(!mpz_cmp_ui(r, 279));

        fp_mul(r, a, b);
        assert(!mpz_cmp_ui(r, 283));

        fp_div(r, a, b);
        assert(!mpz_cmp_ui(r, 11));

        fp_div(r, b, a);
        assert(!mpz_cmp_ui(r, 196));

        fp_clear(r);
    }

    fp_clear(a);
    fp_clear(b);

    assert(!fp_char_clear());
}

int main() {
    test_simple();
    test_modulo();
    return 0;
}