#include <stdio.h>

#include "ec_mont.h"
#include "testing.h"

fp2_t A, C;
fp2_t A24_plus, C24;
point_t P, Q, K, PQd;

void init_test_variables() {
    fp2_init(&A);
    fp2_init(&C);
    fp2_init(&A24_plus);
    fp2_init(&C24);
    point_init(&P);
    point_init(&Q);
    point_init(&K);
    point_init(&PQd);
}

void clear_test_variables() {
    fp2_clear(&A);
    fp2_clear(&C);
    fp2_clear(&A24_plus);
    fp2_clear(&C24);
    point_clear(&P);
    point_clear(&Q);
    point_clear(&K);
    point_clear(&PQd);
}

void set_params_testp431() {
    // Clear just to be sure
    global_fpchar_clear();
    // p + 1 = 432 = 2^4 * 3^3
    global_fpchar_setup_uint(431);

    // Initialize the Elliptic Curve E: y^2 = x^3 + 6x^2 + x
    // (A : C) are projective coordinates of (a : 1) = (6 : 1)
    // We additionally define A24+, C24 and A24- as:
    // A24+: A + 2C
    // C24:  4C
    // A24-: A - 2C
    // This constants are used to speed up the calculations inside 
    // certain functions e.g. xDBL algorithm
    fp2_set_uint(A, 6);
    fp2_set_uint(C, 1);

    // (A, C) = (6, 1)
    // (A24+ : C24) ~ (A + 2C : 4C)
    // This projective pair is the representant of the variable (a + 2)/4
    A24p_from_A(A24_plus, C24, A, C);
}

void set_params_testp139() {
    global_fpchar_clear();
    // p + 1 = 140 = 2^2 * 5 * 7
    global_fpchar_setup_uint(139);

    // E: y^2 = x^3 + x
    // A = 0 => Ax^2 = 0
    fp2_set_uint(A, 6);
    fp2_set_uint(C, 1);

    A24p_from_A(A24_plus, C24, A, C);
}

void test_point_set_is_immutable() {
    fp2_set_str(P->X, "2 + 5*i");
    // Q = P
    point_set(Q, P);
    // Modify Q, but not P
    fp2_set_str(Q->X, "3 + 6*i");

    // Check if point P was modified
    // if yes, then point_set is a shallow copy and not a deepcopy
    CHECK(!mpz_cmp_ui(P->X->a, 2) && !mpz_cmp_ui(P->X->b, 5));
}

void test_xADD_small() {

    // x(P) = XP/1 = 259 + 271i
    point_set_str_x(P, "271*i + 259");
    fp2_print_uint(P->X, "xP");

    // x(Q) = XQ/1 = 262 + 335i
    point_set_str_x(Q, "335*i + 262");
    fp2_print_uint(Q->X, "xQ");

    // x(P - Q) = xPQdiff = 143 + 411i
    point_set_str_x(PQd, "411*i + 143");
    fp2_print_uint(PQd->X, "xP-Q");

    point_t PQsum;
    point_init(&PQsum);

    // Perform the addition:
    xADD(PQsum, P, Q, PQd);

    // x = X(P+Q)/Z(P+Q) = 61 + 184 * i = x(PQ+)
    point_normalize_coords(PQsum);
    fp2_print_uint(PQsum->X, "xP+Q");

    CHECK(!mpz_cmp_ui(PQsum->X->a, 416) && !mpz_cmp_ui(PQsum->X->b, 106));

    point_clear(&PQsum);
}

void test_xDBL_small() {

    // PX = 292 + 15i, zP = 1
    point_set_str_x(P, "15*i + 292");

    fp2_print_uint(P->X, "xP");
    fp2_print_uint(P->Z, "zP");

    // Q = [2]P
    { 
        xDBL(Q, P, A24_plus, C24);

        // Answer in affine coordinates: x(Q) = 61 + 184 * i
        // Which is equal to projective coordinates (this implementation)
        // X(Q) = 157 + 180 * i
        // Z(Q) = 65 + 358 * i
        fp2_print_uint(Q->X, "XQ");
        fp2_print_uint(Q->Z, "ZQ");

        CHECK(!mpz_cmp_ui(Q->X->a, 157) && !mpz_cmp_ui(Q->X->b, 180));
        CHECK(!mpz_cmp_ui(Q->Z->a, 65) && !mpz_cmp_ui(Q->Z->b, 358));

        // X(Q)/Z(Q) = 61 + 184 * i = x(Q)
        point_normalize_coords(Q);
        fp2_print_uint(Q->X, "xQ");
        CHECK(!mpz_cmp_ui(Q->X->a, 61) && !mpz_cmp_ui(Q->X->b, 184));
    }
}

void test_criss_cross_small() {
    fp2_t x, y, z, w, sum, diff;
    fp2_init(&x); fp2_init(&y); fp2_init(&z); fp2_init(&w); fp2_init(&sum); fp2_init(&diff); 

    fp2_set_str(x, "416*i + 175");
    fp2_set_str(y, "112*i + 179");
    fp2_set_str(z, "235*i + 107");
    fp2_set_str(w, "183*i + 197");

    fp2_print_uint(x, "x");
    fp2_print_uint(y, "y");
    fp2_print_uint(z, "z");
    fp2_print_uint(w, "w");
    
    criss_cross(sum, diff, x, y, z, w);

    fp2_print_uint(sum, "sum");
    fp2_print_uint(diff, "diff");

    // ad+bc: 367*i + 314
    CHECK(!mpz_cmp_ui(sum->a, 314) && !mpz_cmp_ui(sum->b, 367));
    // ad-bc: 19*i + 425
    CHECK(!mpz_cmp_ui(diff->a, 425) && !mpz_cmp_ui(diff->b, 19));

    fp2_clear(&x); fp2_init(&y); fp2_init(&z); fp2_init(&w); fp2_init(&sum); fp2_init(&diff); 
}


/*
 * @brief Test whether criss-cross is safe for calling with l or r set to any of (x, y, z, w)
 */
void test_criss_cross_argsafe() {
    fp2_t x, y, z, w;
    fp2_init(&x); fp2_init(&y); fp2_init(&z); fp2_init(&w);

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

    fp2_clear(&x); fp2_init(&y); fp2_init(&z); fp2_init(&w);
}

void test_xLADDER3PT_small() {
    point_set_str_x(P, "271*i + 259");
    point_set_str_x(Q, "335*i + 262");
    point_set_str_x(PQd, "411*i + 143");
    long int n = 87;

    // P contains the result, other points modified
    xLADDER3PT_int(P, Q, PQd, n, A24_plus, C24);
    point_normalize_coords(P);
    
    fp2_print_uint(P->X, "xP");
    CHECK(!mpz_cmp_ui(P->X->a, 360) && !mpz_cmp_ui(P->X->b, 45));

    point_set_str_x(P, "271*i + 259");
    point_set_str_x(Q, "335*i + 262");
    point_set_str_x(PQd, "411*i + 143");

    mpz_t m;
    mpz_init_set_ui(m, n);
    xLADDER3PT(P, Q, PQd, m, A24_plus, C24);
    point_normalize_coords(P);

    fp2_print_uint(P->X, "xP");
    CHECK(!mpz_cmp_ui(P->X->a, 360) && !mpz_cmp_ui(P->X->b, 45));

    mpz_clear(m);
}

void test_point_normalize_coords() {
    fp2_set_str(P->X, "395*i + 201");
    fp2_set_str(P->Z, "272*i + 286");

    // XP = XP / ZP; ZP = 1
    point_normalize_coords(P);
    fp2_print_uint(P->X, "XP");
    fp2_print_uint(P->Z, "ZP");

    CHECK(!mpz_cmp_ui(P->X->a, 95) && !mpz_cmp_ui(P->X->b, 12));
}


void test_KPS() {
    point_set_str_x(K, "101*i + 20");
    fp2_print_uint(K->X, "xK");

    int deg = 7;
    // we do not include last point [deg]K = E(0) which is point at inf.
    const size_t n = KPS_DEG2SIZE(deg);

    point_t kpt[n];
    for (size_t i = 0; i < n; i++) {
        point_init(&kpt[i]);
    }

    KPS(kpt, n, K, A24_plus, C24);

    // Normalize the coordinates for each point
    for (size_t i = 0; i < n; i++) {
        CHECK(!fp2_is_zero(kpt[i]->Z) && "Ki->Z is zero");
        point_normalize_coords(kpt[i]);
    }

    // [ (j+1)*K for j in range(3) ] == [101*i + 20, 82*i + 16, 106*i + 124]
    CHECK(!mpz_cmp_ui(kpt[0]->X->a, 20) && !mpz_cmp_ui(kpt[0]->X->b, 101));
    CHECK(!mpz_cmp_ui(kpt[1]->X->a, 16) && !mpz_cmp_ui(kpt[1]->X->b, 82));
    CHECK(!mpz_cmp_ui(kpt[2]->X->a, 124) && !mpz_cmp_ui(kpt[2]->X->b, 106));

    for (size_t i = 0; i < n; i++) {
        point_clear(&kpt[i]);
    }
}

void test_xISOG_point() {
    // K is the kernel generator of the isogeny
    point_set_str_x(K, "92*i + 97");
    fp2_print_uint(K->X, "xK");

    // P is the point sent through the isogeny
    point_set_str_x(P, "6*i + 114");
    fp2_print_uint(P->X, "xP");

    // isogeny of degree 5 for kernel <K>
    const int degree = 5;
    const size_t n = KPS_DEG2SIZE(degree);
    point_t kpts[n];
    for (size_t i = 0; i < n; i++) point_init(&kpts[i]);
    printf("Size; %lu\n", n);

    KPS(kpts, n, K, A24_plus, C24);
    prepare_kernel_points(kpts, n);

    xISOG_point(Q, kpts, n, P);
    point_normalize_coords(Q);
    fp2_print_uint(Q->X, "xφ(P)");

    for (size_t i = 0; i < n; i++) point_clear(&kpts[i]);
}


int main() {
    init_test_variables();

    // integrity tests
    TEST_RUN(test_point_set_is_immutable());

    // p = 431 tests
    set_params_testp431();

    TEST_RUN(test_point_normalize_coords());
    TEST_RUN(test_xDBL_small());
    TEST_RUN(test_xADD_small());
    TEST_RUN(test_criss_cross_small());
    TEST_RUN(test_criss_cross_argsafe());
    TEST_RUN(test_xLADDER3PT_small());

    // p = 139
    set_params_testp139();
    TEST_RUN(test_KPS());
    TEST_RUN(test_xISOG_point());

    TEST_RUNS_END;

    clear_test_variables();
    return 0;
}