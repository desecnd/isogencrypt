#include <gmp.h>
#include <stdio.h>

#include "fp.h"
#include "isog_mont.h"
#include "fp2.h"
#include "testing.h"

fp2_t A, C;
fp2_t A24p, C24;
point_t P, Q, K, PQd;

void init_test_variables() {
    fp2_init(&A);
    fp2_init(&C);
    fp2_init(&A24p);
    fp2_init(&C24);
    point_init(&P);
    point_init(&Q);
    point_init(&K);
    point_init(&PQd);
}

void clear_test_variables() {
    fp2_clear(&A);
    fp2_clear(&C);
    fp2_clear(&A24p);
    fp2_clear(&C24);
    point_clear(&P);
    point_clear(&Q);
    point_clear(&K);
    point_clear(&PQd);
    fpchar_clear_if_set();
}

void set_params_testp431() {
    // Clear just to be sure
    fpchar_clear_if_set();

    // p + 1 = 432 = 2^4 * 3^3
    assert(0 == fpchar_setup_uint(431));

    // Initialize the Elliptic Curve E: y^2 = x^3 + 6x^2 + x
    // (A : C) are projective coordinates of (a : 1) = (6 : 1)
    // We additionally define A24+, C24:
    // A24+: A + 2C
    // C24:  4C
    // This constants are used to speed up the calculations inside
    // certain functions e.g. xDBL algorithm
    fp2_set_uint(A, 6);
    fp2_set_uint(C, 1);

    // (A, C) = (6, 1)
    // (A24+ : C24) ~ (A + 2C : 4C)
    // This projective pair is the representant of the variable (a + 2)/4
    A24p_from_A(A24p, C24, A, C);
}

void set_params_testp139() {
    // p + 1 = 140 = 2^2 * 5 * 7
    fpchar_clear_if_set();
    assert(0 == fpchar_setup_uint(139));

    // E: y^2 = x^3 + x
    // A = 0 => Ax^2 = 0
    fp2_set_uint(A, 6);
    fp2_set_uint(C, 1);

    A24p_from_A(A24p, C24, A, C);
}

void test_criss_cross_small() {
    fp2_t x, y, z, w, sum, diff;
    fp2_init(&x);
    fp2_init(&y);
    fp2_init(&z);
    fp2_init(&w);
    fp2_init(&sum);
    fp2_init(&diff);

    fp2_set_str(x, "416*i + 175");
    fp2_set_str(y, "112*i + 179");
    fp2_set_str(z, "235*i + 107");
    fp2_set_str(w, "183*i + 197");

    fp2_print(x, "x");
    fp2_print(y, "y");
    fp2_print(z, "z");
    fp2_print(w, "w");

    criss_cross(sum, diff, x, y, z, w);

    fp2_print(sum, "xw+yz");
    fp2_print(diff, "xw-yz");

    // ad+bc: 367*i + 314
    CHECK(!mpz_cmp_ui(sum->a, 314) && !mpz_cmp_ui(sum->b, 367));
    // ad-bc: 19*i + 425
    CHECK(!mpz_cmp_ui(diff->a, 425) && !mpz_cmp_ui(diff->b, 19));

    fp2_clear(&x);
    fp2_clear(&y);
    fp2_clear(&z);
    fp2_clear(&w);
    fp2_clear(&sum);
    fp2_clear(&diff);
}

/*
 * @brief Test whether criss-cross is safe for calling with l or r set to any of
 * (x, y, z, w)
 */
void test_criss_cross_argsafe() {
    fp2_t x, y, z, w;
    fp2_init(&x);
    fp2_init(&y);
    fp2_init(&z);
    fp2_init(&w);

    fp2_set_str(x, "416*i + 175");
    fp2_set_str(y, "112*i + 179");
    fp2_set_str(z, "235*i + 107");
    fp2_set_str(w, "183*i + 197");

    criss_cross(x, y, x, y, z, w);
    CHECK(!mpz_cmp_ui(x->a, 314) && !mpz_cmp_ui(x->b, 367));
    CHECK(!mpz_cmp_ui(y->a, 425) && !mpz_cmp_ui(y->b, 19));

    fp2_set_str(x, "416*i + 175");
    fp2_set_str(y, "112*i + 179");

    criss_cross(z, w, x, y, z, w);
    CHECK(!mpz_cmp_ui(z->a, 314) && !mpz_cmp_ui(z->b, 367));
    CHECK(!mpz_cmp_ui(w->a, 425) && !mpz_cmp_ui(w->b, 19));

    fp2_set_str(z, "235*i + 107");
    fp2_set_str(w, "183*i + 197");

    criss_cross(x, w, x, y, z, w);
    CHECK(!mpz_cmp_ui(x->a, 314) && !mpz_cmp_ui(x->b, 367));
    CHECK(!mpz_cmp_ui(w->a, 425) && !mpz_cmp_ui(w->b, 19));

    fp2_set_str(x, "416*i + 175");
    fp2_set_str(w, "183*i + 197");

    criss_cross(z, y, x, y, z, w);
    CHECK(!mpz_cmp_ui(z->a, 314) && !mpz_cmp_ui(z->b, 367));
    CHECK(!mpz_cmp_ui(y->a, 425) && !mpz_cmp_ui(y->b, 19));

    fp2_clear(&x);
    fp2_clear(&y);
    fp2_clear(&z);
    fp2_clear(&w);
}

void test_ISOG2e() {
    // Point of order 16, does not lay above (0, 0)
    point_set_str_x(K, "33*i + 429");
    fp2_print(K->X, "xK");

    point_set_str_x(P, "158*i + 183");
    fp2_print(P->X, "xP");

    fp2_t A24p_, C24_, A_, C_, a_;
    fp2_init(&A24p_);
    fp2_init(&C24_);
    fp2_init(&a_);
    fp2_init(&A_);
    fp2_init(&C_);

    point_t push_points[2] = {P, NULL};

    // Calculate codomain and push point P through the isogeny of degree 2^4
    // = 16.
    ISOG2e(A24p_, C24_, A24p, C24, K, 4, push_points);

    point_normalize_coords(P);
    fp2_print(P->X, "xφ(P)");

    // xphi(P) = 69*i + 48
    CHECK(!mpz_cmp_ui(P->X->a, 48) && !mpz_cmp_ui(P->X->b, 69));

    A_from_A24p(A_, C_, A24p_, C24_);
    fp2_div_unsafe(a_, A_, C_);
    fp2_print(a_, "aφ(E)");
    // aphi(E) = 201
    CHECK(!mpz_cmp_ui(a_->a, 201) && !mpz_cmp_ui(a_->b, 0));

    fp2_clear(&A24p_);
    fp2_clear(&C24_);
    fp2_clear(&a_);
    fp2_clear(&A_);
    fp2_clear(&C_);
}

// ---------------------
// Testcases for p = 139
// ---------------------

void test_KPS() {
    point_set_str_x(K, "101*i + 20");
    fp2_print(K->X, "xK");

    int deg = 7;
    // we do not include last point [deg]K = E(0) which is point at inf.
    const size_t n = KPS_DEG2SIZE(deg);

    point_t kpt[n];
    for (size_t i = 0; i < n; i++) {
        point_init(&kpt[i]);
    }

    KPS(kpt, n, K, A24p, C24);

    // Normalize the coordinates for each point
    for (size_t i = 0; i < n; i++) {
        CHECK(!fp2_is_zero(kpt[i]->Z) && "Ki->Z is zero");
        point_normalize_coords(kpt[i]);
    }

    // [ (j+1)*K for j in range(3) ] == [101*i + 20, 82*i + 16, 106*i + 124]
    fp2_print(kpt[0]->X, "xK1");
    CHECK(!mpz_cmp_ui(kpt[0]->X->a, 20) && !mpz_cmp_ui(kpt[0]->X->b, 101));
    fp2_print(kpt[1]->X, "xK2");
    CHECK(!mpz_cmp_ui(kpt[1]->X->a, 16) && !mpz_cmp_ui(kpt[1]->X->b, 82));
    fp2_print(kpt[2]->X, "xK3");
    CHECK(!mpz_cmp_ui(kpt[2]->X->a, 124) && !mpz_cmp_ui(kpt[2]->X->b, 106));

    for (size_t i = 0; i < n; i++) {
        point_clear(&kpt[i]);
    }
}

void test_xISOG_and_aISOG() {
    // K is the kernel generator of the isogeny
    point_set_str_x(K, "77*i + 38");
    point_printx(K, "xK");

    // P is the point sent through the isogeny
    point_set_str_x(P, "32*i + 42");
    point_printx(P, "xP");

    // isogeny of degree 5 for kernel <K>
    const int degree = 5;
    const size_t n = KPS_DEG2SIZE(degree);
    point_t kpts[n];
    for (size_t i = 0; i < n; i++)
        point_init(&kpts[i]);
    printf("deg: %d\n", degree);
    printf("n: %lu\n", n);

    KPS(kpts, n, K, A24p, C24);

    fp2_t phiA, phiC, phi_a;
    fp2_init(&phiA);
    fp2_init(&phiC);
    fp2_init(&phi_a);

    // Check codomain Curve value
    aISOG_curve_KPS(phiA, phiC, A24p, C24, kpts, n);
    fp2_div_unsafe(phi_a, phiA, phiC);
    fp2_print(phi_a, "aφ(K)");
    CHECK(!mpz_cmp_ui(phi_a->a, 85) && !mpz_cmp_ui(phi_a->b, 76));

    // Run the same computation, make sure the result is equal
    aISOG_curve(phiA, phiC, A24p, C24, K, degree);
    fp2_div_unsafe(phi_a, phiA, phiC);
    CHECK(!mpz_cmp_ui(phi_a->a, 85) && !mpz_cmp_ui(phi_a->b, 76));

    fp2_clear(&phiA);
    fp2_clear(&phiC);
    fp2_clear(&phi_a);

    prepare_kernel_points(kpts, n);
    xISOG_odd(Q, kpts, n, P);
    point_normalize_coords(Q);
    fp2_print(Q->X, "xφ(P)");

    CHECK(!mpz_cmp_ui(Q->X->a, 46) && !mpz_cmp_ui(Q->X->b, 88));

    for (size_t i = 0; i < n; i++)
        point_clear(&kpts[i]);
}

// TODO: reform the test
void test_ISOG_chain_odd() {

    // Point of order 35 on the curve E
    point_set_str_x(K, "108*i + 136");
    fp2_print(K->X, "xK");

    fp2_t A_, C_, a_;
    fp2_init(&A_);
    fp2_init(&C_), fp2_init(&a_);

    pprod_t deg;
    pprod_init(&deg);

    unsigned int primes[2] = {5, 7};
    pprod_set_array(deg, primes, 2);

    point_t push_points[] = {NULL, NULL};
    ISOG_chain(A_, C_, A24p, C24, K, deg, push_points);

    // aφ(K): 102*i + 73
    A_from_A24p(A_, C_, A_, C_);
    fp2_div_unsafe(a_, A_, C_);
    fp2_print(a_, "aφ(K)");
    CHECK(!mpz_cmp_ui(a_->a, 73) && !mpz_cmp_ui(a_->b, 102));

    fp2_clear(&A_);
    fp2_clear(&C_);
    fp2_clear(&a_);
    pprod_clear(&deg);
}

void test_xISOG2_and_aISOG2() {
    // Point of order 2, x != 0
    point_set_str_x(K, "100*i + 136");
    point_printx(K, "xK2");

    // Point of order 140
    point_set_str_x(P, "70*i + 36");
    point_printx(P, "xP");

    fp2_t A_, C_, a_;
    fp2_init(&A_);
    fp2_init(&C_), fp2_init(&a_);

    // Calculate codomain of the 2-isogeny curve
    // a' = 37*i + 73
    aISOG2(A_, C_, K);
    fp2_div_unsafe(a_, A_, C_);
    fp2_print(a_, "aE2");
    CHECK(!mpz_cmp_ui(a_->a, 73) && !mpz_cmp_ui(a_->b, 37));

    // Calculate codomain of the 2-isogeny curve (in xDBL form)
    // a' = 44*i + 123
    aISOG2_24p(A_, C_, K);
    fp2_div_unsafe(a_, A_, C_);
    fp2_print(a_, "aE2(24p)");
    CHECK(!mpz_cmp_ui(a_->a, 123) && !mpz_cmp_ui(a_->b, 44));

    // Calculate the x-coordinate of the image of P under the 2-isogeny
    // xP' = 128*i
    xISOG2_unsafe(Q, K, P);
    point_normalize_coords(Q);
    point_printx(Q, "xφ(P)");
    CHECK(!mpz_cmp_ui(Q->X->a, 0) && !mpz_cmp_ui(Q->X->b, 128));

    // Calculate the x-coordinate using prepared Kernel
    prepare_isog2_kernel(K);
    xISOG2_prep(Q, K, P);
    point_normalize_coords(Q);
    CHECK(!mpz_cmp_ui(Q->X->a, 0) && !mpz_cmp_ui(Q->X->b, 128));

    fp2_clear(&A_);
    fp2_clear(&C_), fp2_clear(&a_);
}

void test_ISOG_chain() {
    // Point of order 140, does not lay over (0, 0)
    point_set_str_x(K, "34*i + 99");
    point_printx(K, "xK");

    point_set_str_x(P, "8*i + 137");
    point_printx(P, "xP");

    fp2_t A_, C_, a_, A24p_, C24_;
    fp2_init(&A24p_);
    fp2_init(&C24_);
    fp2_init(&A_);
    fp2_init(&C_), fp2_init(&a_);

    pprod_t deg;
    pprod_init(&deg);

    unsigned int factors[] = {4, 5, 7};
    pprod_set_array(deg, factors, 3);

    // Two more NULLs for usage by ISOG_chain
    point_t push_points[3] = {P, NULL, NULL};

    // Return in 24p form
    ISOG_chain(A24p_, C24_, A24p, C24, K, deg, push_points);
    A_from_A24p(A_, C_, A24p_, C24_);
    fp2_div_unsafe(a_, A_, C_);

    fp2_print(a_, "aφ(K)");
    CHECK(!mpz_cmp_ui(a_->a, 0) && !mpz_cmp_ui(a_->b, 49));

    point_normalize_coords(P);
    point_printx(P, "xφ(P)");
    CHECK(!mpz_cmp_ui(P->X->a, 68) && !mpz_cmp_ui(P->X->b, 114));

    fp2_clear(&A24p_);
    fp2_clear(&C24_);
    fp2_clear(&A_);
    fp2_clear(&C_), fp2_clear(&a_);
    pprod_clear(&deg);
}

int main() {
    init_test_variables();

    // p = 431 tests
    set_params_testp431();

    TEST_RUN_SILENT(test_criss_cross_argsafe());
    TEST_RUN(test_criss_cross_small());
    TEST_RUN(test_ISOG2e());

    // p = 139 tests
    set_params_testp139();

    TEST_RUN(test_KPS());
    TEST_RUN(test_xISOG_and_aISOG());
    TEST_RUN(test_ISOG_chain_odd());
    TEST_RUN(test_xISOG2_and_aISOG2());
    TEST_RUN(test_ISOG_chain());

    clear_test_variables();

    TEST_RUNS_END;
}