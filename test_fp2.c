#include <stdio.h>

#include "fp2.h"
#include "testing.h"

void test_failed_once() {
    CHECK(!global_fpchar_setup_uint(431));

    fp2_t x, y;
    fp2_init(&x);
    fp2_init(&y);

    CHECK(fp2_is_zero(x) && fp2_is_zero(y));

    // x = 291 + 15*i
    fp_set_uint(x->a, 291);
    fp_set_uint(x->b, 15);
    fp2_print_uint(x, "x");

    // y = 293 + 15*i
    fp_set_uint(y->a, 293);
    fp_set_uint(y->b, 15);
    fp2_print_uint(y, "y");

    // Square-ing test that failed
    {
        fp2_t r; 
        fp2_init(&r);

        fp2_sq_unsafe(r, x);
        fp2_print_uint(r, "x^2");
        CHECK(!mpz_cmp_ui(r->a, 411));
        CHECK(!mpz_cmp_ui(r->b, 110));

        fp2_sq_unsafe(r, y);
        fp2_print_uint(r, "y^2");
        CHECK(!mpz_cmp_ui(r->a, 286));
        CHECK(!mpz_cmp_ui(r->b, 170));

        fp2_clear(&r);
    }

    CHECK(!global_fpchar_clear());
}

void test_basic_arithmetic() {
    printf("[*] Running Fp^2 arithmetic test for char p = 431...\n");

    // init characteristic
    CHECK(!global_fpchar_setup_uint(431));

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
    CHECK(!mpz_cmp_ui(r->a, 68));
    CHECK(!mpz_cmp_ui(r->b, 84));

    // r <- x - y
    // r = (27 + 59i) - (47 + 25i)
    // r = (-26 + 34i) 
    // r = 405 + 34i (mod 431)
    fp2_sub(r, x, y);
    printf("x - y: %ld, %ldi\n", mpz_get_ui(r->a), mpz_get_ui(r->b));
    CHECK(!mpz_cmp_ui(r->a, 405));
    CHECK(!mpz_cmp_ui(r->b, 34));

    // r <- x * y
    // r = (27 + 59i) * (47 + 25i)
    // r = (-488 + 3298i)
    // r = 374 + 281i (mod 431)
    fp2_mul_unsafe(r, x, y);
    printf("x * y: %ld, %ldi\n", mpz_get_ui(r->a), mpz_get_ui(r->b));
    CHECK(!mpz_cmp_ui(r->a, 374));
    CHECK(!mpz_cmp_ui(r->b, 281));

    // r <- x^2
    // r = (27 + 59i) * (27 + 59i)
    // r = (-3040 + 2478)
    // r = 408 + 323i (mod 431)
    fp2_sq_unsafe(r, x);
    printf("x ^ 2: %ld, %ldi\n", mpz_get_ui(r->a), mpz_get_ui(r->b));
    CHECK(!mpz_cmp_ui(r->a, 408));
    CHECK(!mpz_cmp_ui(r->b, 323));

    fp2_clear(&x);
    fp2_clear(&y);
    fp2_clear(&r);

    CHECK(!global_fpchar_clear());
}

void test_safe_unsafe_mul_sq() {

    // init characteristic
    CHECK(!global_fpchar_setup_uint(431));

    fp2_t r, x;
    fp2_init(&r); fp2_init(&x);

    // x = 27 + 59i
    fp_set_uint(x->a, 21);
    fp_set_uint(x->b, 59);
    fp2_print_uint(x, "x");

    fp2_set(r, x);
    fp2_print_uint(r, "r");

    fp2_mul_unsafe(r, x, x);
    fp2_print_uint(r, "(unsafe) r := x * x");
    CHECK(!mpz_cmp_ui(r->a, 408) && !mpz_cmp_ui(r->b, 323));

    fp2_set(r, x);
    fp2_mul_safe(r, x);
    fp2_print_uint(r, "(safe) r *= x");
    CHECK(!mpz_cmp_ui(r->a, 408) && !mpz_cmp_ui(r->b, 323));

    fp2_set(r, x);
    fp2_mul_safe(r,  r);
    fp2_print_uint(r, "(safe) r *= r");
    CHECK(!mpz_cmp_ui(r->a, 408) && !mpz_cmp_ui(r->b, 323));

    fp2_set(r, x);
    fp2_sq_unsafe(r, x);
    fp2_print_uint(r, "(unsafe-sq) r := x^2");
    CHECK(!mpz_cmp_ui(r->a, 408) && !mpz_cmp_ui(r->b, 323));

    fp2_set(r, x);
    fp2_sq_safe(r);
    fp2_print_uint(r, "(safe-sq) r ^= 2");
    CHECK(!mpz_cmp_ui(r->a, 408) && !mpz_cmp_ui(r->b, 323));

    // Make sure that unsafe is wrong when called wrong?
    // Currently mul contains CHECKion that prevents this check
    // fp2_set(r, x);
    // fp2_mul_unsafe(r, r, r);
    // fp2_print_uint(r, "(unsafe) r := r * r");
    // CHECK(mpz_cmp_ui(r->a, 408) || mpz_cmp_ui(r->b, 323));

    fp2_clear(&x);
    fp2_clear(&r);

    CHECK(!global_fpchar_clear());
}

void test_inv_div() {
    CHECK(!global_fpchar_setup_uint(431));

    fp2_t x, r;
    fp2_init(&x);
    fp2_init(&r);

    {
        // x = 27 + 59i
        fp_set_uint(x->a, 21);
        fp_set_uint(x->b, 59);
        fp2_print_uint(x, "x");

        // r = ~x = 221 + 159i
        fp2_inv_unsafe(r, x);
        fp2_print_uint(r, "~x");
        CHECK(!mpz_cmp_ui(r->a, 221) && !mpz_cmp_ui(r->b, 159));

        // r * x == 1
        fp2_mul_safe(r, x);
        fp2_print_uint(r, "~x * x");
        CHECK(!mpz_cmp_ui(r->a, 1) && !mpz_cmp_ui(r->b, 0));

    }

    {
        // x = 27 + 59i
        fp_set_uint(x->a, 270);
        fp_set_uint(x->b, 240);
        fp2_print_uint(x, "x");

        // r = ~x = 221 + 159i
        fp2_inv_unsafe(r, x);
        fp2_print_uint(r, "~x");
        CHECK(!mpz_cmp_ui(r->a, 11) && !mpz_cmp_ui(r->b, 86));

        // r * x == 1
        fp2_mul_safe(r, x);
        fp2_print_uint(r, "~x * x");
        CHECK(!mpz_cmp_ui(r->a, 1) && !mpz_cmp_ui(r->b, 0));

    }

    fp2_clear(&x);
    fp2_clear(&r);

    CHECK(!global_fpchar_clear());
}

void test_set_str() {
    CHECK(!global_fpchar_setup_uint(431));

    fp2_t x;
    fp2_init(&x);

    fp2_set_str(x, "416*i + 175");
    fp2_print_uint(x, "x: ");
    CHECK(!mpz_cmp_ui(x->a, 175) && !mpz_cmp_ui(x->b, 416));

    fp2_set_str(x, "416 + 175 * i");
    CHECK(!mpz_cmp_ui(x->a, 416) && !mpz_cmp_ui(x->b, 175));

    fp2_set_str(x, "0xff + 23 * i");
    CHECK(!mpz_cmp_ui(x->a, 255) && !mpz_cmp_ui(x->b, 23));

    CHECK(0 != fp2_set_str(x, "10"));
    CHECK(0 != fp2_set_str(x, "10 * i"));

    fp2_clear(&x);
    CHECK(!global_fpchar_clear());
}

void test_fp2_mul_int() {
    CHECK(!global_fpchar_setup_uint(431));

    fp2_t r, x;
    fp2_init(&r); fp2_init(&x);
    fp2_set_str(x, "224*i + 192");
    fp2_print_uint(x, "x");

    // [224i + 192] * 365 == 301*i + 258 (mod 431)
    fp2_mul_int(r, x, 365);
    fp2_print_uint(r, "365x");
    CHECK(!mpz_cmp_ui(r->a, 258) && !mpz_cmp_ui(r->b, 301));

    // [224i + 192] * -5 == 173*i + 333 (mod 431)
    fp2_mul_int(r, x, -5);
    fp2_print_uint(r, "-5x");
    CHECK(!mpz_cmp_ui(r->a, 333) && !mpz_cmp_ui(r->b, 173));

    // [145*i + 309] * 1234 == 145*i + 309  (mod 431)
    fp2_mul_int(r, x, 1234);
    fp2_print_uint(r, "1234x");
    CHECK(!mpz_cmp_ui(r->a, 309) && !mpz_cmp_ui(r->b, 145));

    fp2_clear(&r); fp2_clear(&x);
    CHECK(!global_fpchar_clear());
}

int main() {

    TEST_RUN(test_set_str());
    TEST_RUN(test_failed_once());
    TEST_RUN(test_basic_arithmetic());
    TEST_RUN(test_safe_unsafe_mul_sq());
    TEST_RUN(test_inv_div());
    TEST_RUN(test_fp2_mul_int());
    TEST_RUNS_END;

    return 0;
}