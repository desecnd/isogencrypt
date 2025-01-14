#include <stdio.h>
#include <assert.h>

#include "testing.h"
#include "fp2.h"

void test_small_arithmetic() {
    CHECK(!global_fpchar_setup_uint(431));

    fp_t r;
    fp_init(r);

    CHECK(fp_is_zero(r));

    fp_set_uint(r, 16); // r == 16
    CHECK(!mpz_cmp_ui(r, 16));

    fp_sqrt(r, r); // r == 4
    // TODO: incosistent api, change later
    CHECK(!mpz_cmp_ui(r, 4));

    fp_add_uint(r, r, 10); // r == 14
    CHECK(!mpz_cmp_ui(r, 14));

    fp_sub_uint(r, r, 15); // r == -1 == 430
    CHECK(!mpz_cmp_ui(r, 430));

    fp_inv(r, r); // r == r ^-1 == 430
    CHECK(!mpz_cmp_ui(r, 430));

    fp_clear(r);

    CHECK(!global_fpchar_clear());
}

void test_modulo_arithmetic() {
    CHECK(!global_fpchar_setup_uint(431));

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
        CHECK(!mpz_cmp_ui(r, 10));

        fp_sub(r, a, b);
        CHECK(!mpz_cmp_ui(r, 152));

        fp_sub(r, b, a);
        CHECK(!mpz_cmp_ui(r, 279));

        fp_mul(r, a, b);
        CHECK(!mpz_cmp_ui(r, 283));

        fp_div(r, a, b);
        CHECK(!mpz_cmp_ui(r, 11));

        fp_div(r, b, a);
        CHECK(!mpz_cmp_ui(r, 196));

        fp_clear(r);
    }

    fp_clear(a);
    fp_clear(b);

    CHECK(!global_fpchar_clear());
}

int main() {

    TEST_RUN(test_small_arithmetic());
    TEST_RUN(test_modulo_arithmetic());
    TEST_RUNS_END;

    return 0;
}