#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fp2.h"
#include "testing.h"

void test_failed_once() {
    CHECK(!fpchar_setup_uint(431));

    fp2_t x, y;
    fp2_init(&x);
    fp2_init(&y);

    CHECK(fp2_is_zero(x) && fp2_is_zero(y));

    // x = 291 + 15*i
    fp_set_uint(x->a, 291);
    fp_set_uint(x->b, 15);
    fp2_print(x, "x");

    // y = 293 + 15*i
    fp_set_uint(y->a, 293);
    fp_set_uint(y->b, 15);
    fp2_print(y, "y");

    // Square-ing test that failed
    {
        fp2_t r; 
        fp2_init(&r);

        fp2_sq_unsafe(r, x);
        fp2_print(r, "x^2");
        CHECK(!mpz_cmp_ui(r->a, 411));
        CHECK(!mpz_cmp_ui(r->b, 110));

        fp2_sq_unsafe(r, y);
        fp2_print(r, "y^2");
        CHECK(!mpz_cmp_ui(r->a, 286));
        CHECK(!mpz_cmp_ui(r->b, 170));

        fp2_clear(&r);
    }

    CHECK(!fpchar_clear());
}

void test_basic_arithmetic() {
    // init characteristic
    CHECK(!fpchar_setup_uint(431));

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

    CHECK(!fpchar_clear());
}

void test_safe_unsafe_mul_sq() {

    // init characteristic
    CHECK(!fpchar_setup_uint(431));

    fp2_t r, x;
    fp2_init(&r); fp2_init(&x);

    // x = 27 + 59i
    fp_set_uint(x->a, 21);
    fp_set_uint(x->b, 59);
    fp2_print(x, "x");

    fp2_set(r, x);
    fp2_print(r, "r");

    fp2_mul_unsafe(r, x, x);
    fp2_print(r, "(unsafe) r := x * x");
    CHECK(!mpz_cmp_ui(r->a, 408) && !mpz_cmp_ui(r->b, 323));

    fp2_set(r, x);
    fp2_mul_safe(r, x);
    fp2_print(r, "(safe) r *= x");
    CHECK(!mpz_cmp_ui(r->a, 408) && !mpz_cmp_ui(r->b, 323));

    fp2_set(r, x);
    fp2_mul_safe(r,  r);
    fp2_print(r, "(safe) r *= r");
    CHECK(!mpz_cmp_ui(r->a, 408) && !mpz_cmp_ui(r->b, 323));

    fp2_set(r, x);
    fp2_sq_unsafe(r, x);
    fp2_print(r, "(unsafe-sq) r := x^2");
    CHECK(!mpz_cmp_ui(r->a, 408) && !mpz_cmp_ui(r->b, 323));

    fp2_set(r, x);
    fp2_sq_safe(r);
    fp2_print(r, "(safe-sq) r ^= 2");
    CHECK(!mpz_cmp_ui(r->a, 408) && !mpz_cmp_ui(r->b, 323));

    // Make sure that unsafe is wrong when called wrong?
    // Currently mul contains CHECKion that prevents this check
    // fp2_set(r, x);
    // fp2_mul_unsafe(r, r, r);
    // fp2_print_uint(r, "(unsafe) r := r * r");
    // CHECK(mpz_cmp_ui(r->a, 408) || mpz_cmp_ui(r->b, 323));

    fp2_clear(&x);
    fp2_clear(&r);

    CHECK(!fpchar_clear());
}

void test_inv_div() {
    CHECK(!fpchar_setup_uint(431));

    fp2_t x, r;
    fp2_init(&x);
    fp2_init(&r);

    {
        // x = 27 + 59i
        fp_set_uint(x->a, 21);
        fp_set_uint(x->b, 59);
        fp2_print(x, "x");

        // r = ~x = 221 + 159i
        fp2_inv_unsafe(r, x);
        fp2_print(r, "~x");
        CHECK(!mpz_cmp_ui(r->a, 221) && !mpz_cmp_ui(r->b, 159));

        // r * x == 1
        fp2_mul_safe(r, x);
        fp2_print(r, "~x * x");
        CHECK(!mpz_cmp_ui(r->a, 1) && !mpz_cmp_ui(r->b, 0));

    }

    {
        // x = 27 + 59i
        fp_set_uint(x->a, 270);
        fp_set_uint(x->b, 240);
        fp2_print(x, "x");

        // r = ~x = 221 + 159i
        fp2_inv_unsafe(r, x);
        fp2_print(r, "~x");
        CHECK(!mpz_cmp_ui(r->a, 11) && !mpz_cmp_ui(r->b, 86));

        // r * x == 1
        fp2_mul_safe(r, x);
        fp2_print(r, "~x * x");
        CHECK(!mpz_cmp_ui(r->a, 1) && !mpz_cmp_ui(r->b, 0));

    }

    fp2_clear(&x);
    fp2_clear(&r);

    CHECK(!fpchar_clear());
}

void test_set_str() {
    CHECK(!fpchar_setup_uint(431));

    fp2_t x;
    fp2_init(&x);

    // Valid
    CHECK(0 == fp2_set_str(x, "416*i + 175"));
    CHECK(!mpz_cmp_ui(x->a, 175) && !mpz_cmp_ui(x->b, 416));

    CHECK(0 == fp2_set_str(x, "416 + 175 *i"));
    CHECK(!mpz_cmp_ui(x->a, 416) && !mpz_cmp_ui(x->b, 175));

    CHECK(0 == fp2_set_str(x, "0xff + 23 *i"));
    CHECK(!mpz_cmp_ui(x->a, 255) && !mpz_cmp_ui(x->b, 23));

    CHECK(0 == fp2_set_str(x, "0xff"));
    CHECK(!mpz_cmp_ui(x->a, 255) && !mpz_cmp_ui(x->b, 0));

    CHECK(0 == fp2_set_str(x, "0xff*i"));
    CHECK(!mpz_cmp_ui(x->a, 0) && !mpz_cmp_ui(x->b, 255));

    CHECK(0 == fp2_set_str(x, "420      *i"));
    CHECK(!mpz_cmp_ui(x->a, 0) && !mpz_cmp_ui(x->b, 420));

    CHECK(0 == fp2_set_str(x, "7+5*i"));
    CHECK(!mpz_cmp_ui(x->a, 7) && !mpz_cmp_ui(x->b, 5));

    // Invalid
    CHECK(0 != fp2_set_str(x, "i*420"));
    CHECK(0 != fp2_set_str(x, "10 * i + 10 * i"));
    CHECK(0 != fp2_set_str(x, "10 * i"));
    CHECK(0 != fp2_set_str(x, "1000; * i"));

    fp2_clear(&x);
    CHECK(!fpchar_clear());
}

void test_write() {
    CHECK(!fpchar_setup_uint(431));
    fp2_t x;
    fp2_init(&x);

    CHECK(0 == fp2_set_str(x, "416*i + 175"));
    size_t buff_size = fp2_write_size(x);
    CHECK(3 + 5 + 3 + 1 == buff_size);

    char * buffer = (char*) malloc(buff_size);
    fp2_write(x, buffer);
    CHECK(3 + 5 + 3 == strlen(buffer));
    CHECK(0 == strcmp(buffer, "416*i + 175"));

    fp2_t y;
    fp2_init(&y);
    fp2_set_str(y, buffer);
    CHECK(fp2_equal(x, y));

    fp2_clear(&y);
    fp2_clear(&x);
    CHECK(!fpchar_clear());
}

void test_fp2_mul_int() {
    CHECK(!fpchar_setup_uint(431));

    fp2_t r, x;
    fp2_init(&r); fp2_init(&x);
    fp2_set_str(x, "224*i + 192");
    fp2_print(x, "x");

    // [224i + 192] * 365 == 301*i + 258 (mod 431)
    fp2_mul_int(r, x, 365);
    fp2_print(r, "365x");
    CHECK(!mpz_cmp_ui(r->a, 258) && !mpz_cmp_ui(r->b, 301));

    // [224i + 192] * -5 == 173*i + 333 (mod 431)
    fp2_mul_int(r, x, -5);
    fp2_print(r, "-5x");
    CHECK(!mpz_cmp_ui(r->a, 333) && !mpz_cmp_ui(r->b, 173));

    // [145*i + 309] * 1234 == 145*i + 309  (mod 431)
    fp2_mul_int(r, x, 1234);
    fp2_print(r, "1234x");
    CHECK(!mpz_cmp_ui(r->a, 309) && !mpz_cmp_ui(r->b, 145));

    fp2_clear(&r); fp2_clear(&x);
    CHECK(!fpchar_clear());
}

void test_large_numbers() {
    mpz_t p;
    mpz_init(p);
    mpz_set_str(p, "14475d5aeccf245fce0e61716bd33537235ad8c4a76a401a4eb1a0fb9cb477dfb", 16);

    fpchar_setup(p);

    fp2_t x, y, r;
    fp2_init(&x);
    fp2_init(&y);
    fp2_init(&r);

    fp2_set_str(x, "94579799687746926965252261998050766088232049658745702162469909894445951219281*i + 39836810185194070024401365631938781071356172459462816998642032429396657654496");
    fp2_set_str(y, "71596955964097321089638264126073859217891501550880385776393603069176960285685*i + 6005039997071591805159957858046559317574936997882312151258378710713185921308");

    // r = x + y
    fp2_add(r, x, y);
    CHECK(fp2_equal_str(r, "19419196819700907788901412019741448840917441203891366685074144967702566597707*i + 45841850182265661829561323489985340388931109457345129149900411140109843575804"));

    // r = x - y
    fp2_sub(r, x, y);
    CHECK(fp2_equal_str(r, "22982843723649605875613997871976906870340548107865316386076306825268990933596*i + 33831770188122478219241407773892221753781235461580504847383653718683471733188"));

    // r = x * y
    fp2_mul_unsafe(r, x, y);
    CHECK(fp2_equal_str(r, "143643133666762529271905768435498468336623797700446365668130846446718446017746*i + 87376311329126807534668591296448536610706136584366568663296846757602776963064"));

    // r = x; r *= y
    fp2_set(r, x);
    fp2_mul_safe(r, y);
    CHECK(fp2_equal_str(r, "143643133666762529271905768435498468336623797700446365668130846446718446017746*i + 87376311329126807534668591296448536610706136584366568663296846757602776963064"));

    // r = x / y = x * y^-1
    fp2_div_unsafe(r, x, y); 
    CHECK(fp2_equal_str(r, "65466103758765061295109857012386042169348498301553700172743066711782897011018*i + 14327469402166493311065643408383417673951894623189817384874372125280170171495"));

    // r = x^2
    fp2_sq_unsafe(r, x);
    CHECK(fp2_equal_str(r, "97160454713256180032332857022443502045185781961919588929542931057357925067882*i + 122888906885422413092835230938250002121611842268753523538473182063181056580234"));

    // r = x; r **= 2
    fp2_set(r, x);
    fp2_sq_safe(r);
    CHECK(fp2_equal_str(r, "97160454713256180032332857022443502045185781961919588929542931057357925067882*i + 122888906885422413092835230938250002121611842268753523538473182063181056580234"));

    fp2_clear(&x);
    fp2_clear(&y);
    fp2_clear(&r);

    mpz_clear(p);

    fpchar_clear();
}

int main() {
    TEST_RUN(test_set_str());
    TEST_RUN_SILENT(test_write());
    TEST_RUN(test_failed_once());
    TEST_RUN(test_basic_arithmetic());
    TEST_RUN(test_safe_unsafe_mul_sq());
    TEST_RUN(test_inv_div());
    TEST_RUN(test_fp2_mul_int());
    TEST_RUN(test_large_numbers());
    TEST_RUNS_END;
}