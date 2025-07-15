#include <gmp.h>
#include <stdio.h>

#include "ec_mont.h"
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

/*
 * @brief Make sure that calling A24p_from_A and A_from_A24p is safe with 'same'
 * arguments: func(A, C, A, C)
 */
void test_A_A24p_conversion_in_place() {
    fp2_t At, Ct, at;
    fp2_init(&At);
    fp2_init(&Ct);
    fp2_init(&at);

    // p + 1 = 432 = 2^4 * 3^3
    fpchar_setup_uint(431);

    int a = 410;
    fp2_set_uint(At, a);
    fp2_set_uint(Ct, 1);

    // Multiply by random scalar
    fp2_mul_int(At, At, 12341234);
    fp2_mul_int(Ct, Ct, 12341234);

    fp2_div_unsafe(at, At, Ct);
    CHECK(fp2_equal_uint(at, a));

    int a24p = 103;
    CHECK(a24p * 4 - 2 == a);

    // Convert to (a + 2)/4 form
    A24p_from_A(At, Ct, At, Ct);
    fp2_div_unsafe(at, At, Ct);
    CHECK(fp2_equal_uint(at, a24p));

    fp2_mul_int(At, At, 98769876);
    fp2_mul_int(Ct, Ct, 98769876);

    // Convert from (a + 2)/4 form
    A_from_A24p(At, Ct, At, Ct);
    fp2_div_unsafe(at, At, Ct);
    CHECK(fp2_equal_uint(at, a));

    fpchar_clear();
    fp2_clear(&At);
    fp2_clear(&Ct);
    fp2_clear(&at);
}

void test_A_A24p_conversion() {
    // Clear just to be sure
    // p + 1 = 432 = 2^4 * 3^3
    fpchar_setup_uint(431);

    fp2_t a;
    fp2_init(&a);

    // (A: C) = (a : 1) = 6
    fp2_set_uint(A, 6);
    fp2_set_uint(C, 1);

    // Multiply by random scalar
    fp2_mul_int(A, A, 1000);
    fp2_mul_int(C, C, 1000);

    A24p_from_A(A24p, C24, A, C);

    // normalize to affine: a = A24p / C24
    fp2_div_unsafe(a, A24p, C24);

    // (A24p : C24) = (a + 2 : 4) = 2
    CHECK(fp2_equal_uint(a, 2));

    // Retrieve back (A : C) given (A24p : C24)
    A_from_A24p(A, C, A24p, C24);
    fp2_div_unsafe(a, A, C);
    CHECK(fp2_equal_uint(a, 6));

    fp2_clear(&a);
    fpchar_clear();
}

void test_point_set_is_immutable() {
    point_set_str_x(P, "5*i + 2");
    // Q = P
    point_set(Q, P);
    // Modify Q, but not P
    point_set_str_x(Q, "6*i + 3");

    // Check if point P was modified
    // if yes, then point_set is a shallow copy and not a deepcopy
    CHECK(!mpz_cmp_ui(P->X->a, 2) && !mpz_cmp_ui(P->X->b, 5));
}

// ------------------------------
// | SIDH-like prime tests p431
// ------------------------------

void test_point_normalize_coords() {
    fp2_set_str(P->X, "395*i + 201");
    fp2_set_str(P->Z, "272*i + 286");

    fp2_print(P->X, "X");
    fp2_print(P->Z, "Y");

    CHECK(!point_is_normalized(P));

    // XP = XP / ZP; ZP = 1
    point_normalize_coords(P);
    CHECK(point_is_normalized(P));

    fp2_print(P->X, "x");
    CHECK(fp2_equal_str(P->X, "12*i + 95"));
}

void test_xDBL_small() {

    // PX = 292 + 15i, zP = 1
    point_set_str_x(P, "15*i + 292");
    point_printx(P, "xP");

    // Q = [2]P
    xDBL(Q, P, A24p, C24);

    // Answer in affine coordinates: x(Q) = 61 + 184 * i
    // Which is equal to projective coordinates (in this implementation)
    // X(Q) = 157 + 180 * i
    // Z(Q) = 65 + 358 * i
    CHECK(fp2_equal_str(Q->X, "180*i + 157"));
    CHECK(fp2_equal_str(Q->Z, "358*i + 65"));

    // X(Q)/Z(Q) = 61 + 184 * i = x(Q)
    point_normalize_coords(Q);
    fp2_print(Q->X, "x[2]P");

    CHECK(fp2_equal_str(Q->X, "184*i + 61"));
}

void test_xDBLe() {
    point_set_str_x(P, "387*i + 387");
    point_printx(P, "xP");

    // x[2]P: 400*i + 311
    xDBLe(Q, P, A24p, C24, 1);
    point_normalize_coords(Q);

    point_printx(Q, "x[2]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 311) && !mpz_cmp_ui(Q->X->b, 400));

    // x[4]P: 13*i + 67
    xDBLe(Q, P, A24p, C24, 2);
    point_normalize_coords(Q);

    point_printx(Q, "x[4]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 67) && !mpz_cmp_ui(Q->X->b, 13));

    // x[8]P: 213*i + 105
    xDBLe(Q, P, A24p, C24, 3);
    point_normalize_coords(Q);

    point_printx(Q, "x[8]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 105) && !mpz_cmp_ui(Q->X->b, 213));

    // x[2^1235]P: 304*i + 223
    xDBLe(Q, P, A24p, C24, 12345);
    point_normalize_coords(Q);

    point_printx(Q, "x[2^12345]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 223) && !mpz_cmp_ui(Q->X->b, 304));
}

void test_xADD_small() {

    // x(P) = XP/1 = 259 + 271i
    point_set_str_x(P, "271*i + 259");
    fp2_print(P->X, "xP");

    // x(Q) = XQ/1 = 262 + 335i
    point_set_str_x(Q, "335*i + 262");
    fp2_print(Q->X, "xQ");

    // x(P - Q) = xPQdiff = 143 + 411i
    point_set_str_x(PQd, "411*i + 143");
    fp2_print(PQd->X, "xP-Q");

    point_t PQsum;
    point_init(&PQsum);

    // Perform the addition:
    xADD(PQsum, P, Q, PQd);

    // x = X(P+Q)/Z(P+Q) = 61 + 184 * i = x(PQ+)
    point_normalize_coords(PQsum);
    fp2_print(PQsum->X, "xP+Q");

    CHECK(!mpz_cmp_ui(PQsum->X->a, 416) && !mpz_cmp_ui(PQsum->X->b, 106));

    point_clear(&PQsum);
}

void test_xLADDER3PT() {
    point_set_str_x(P, "271*i + 259");
    point_set_str_x(Q, "335*i + 262");
    point_set_str_x(PQd, "411*i + 143");
    long int n = 87;

    // Print the init variables
    point_printx(P, "xP");
    point_printx(Q, "xQ");
    point_printx(PQd, "xP-Q");
    printf("n: %ld\n", n);

    // P contains the result, other points modified
    xLADDER3PT_int(P, Q, PQd, n, A24p, C24);
    point_normalize_coords(P);

    point_printx(P, "x(P+nQ)");
    CHECK(!mpz_cmp_ui(P->X->a, 360) && !mpz_cmp_ui(P->X->b, 45));

    point_set_str_x(P, "271*i + 259");
    point_set_str_x(Q, "335*i + 262");
    point_set_str_x(PQd, "411*i + 143");

    // Make sure that result is the same for not general xLADDER3PT
    mpz_t m;
    mpz_init_set_ui(m, n);
    xLADDER3PT(P, Q, PQd, m, A24p, C24);

    point_normalize_coords(P);
    CHECK(!mpz_cmp_ui(P->X->a, 360) && !mpz_cmp_ui(P->X->b, 45));

    mpz_clear(m);
}

// ---------------------
// Testcases for p = 139
// ---------------------


void test_xLADDER_int() {

    point_set_str_x(P, "108*i + 136");
    point_printx(P, "xP");

    // x[1]P = 108*i + 136
    xLADDER_int(Q, P, 1, A24p, C24);
    point_normalize_coords(Q);

    point_printx(Q, "x[1]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 136) && !mpz_cmp_ui(Q->X->b, 108));

    // x[2]P = 113*i + 131
    xLADDER_int(Q, P, 2, A24p, C24);
    point_normalize_coords(Q);

    point_printx(Q, "x[2]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 131) && !mpz_cmp_ui(Q->X->b, 113));

    // x[3]P = 42*i + 83
    xLADDER_int(Q, P, 3, A24p, C24);
    point_normalize_coords(Q);

    point_printx(Q, "x[3]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 83) && !mpz_cmp_ui(Q->X->b, 42));

    // x[4]P = 47*i + 107
    xLADDER_int(Q, P, 4, A24p, C24);
    point_normalize_coords(Q);

    point_printx(Q, "x[4]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 107) && !mpz_cmp_ui(Q->X->b, 47));
}

void test_xLADDER() {
    point_set_str_x(P, "7*i + 97");
    fp2_print(P->X, "xP");

    mpz_t m;
    // 1. m = 2^80
    mpz_init_set_str(m, "100000000000000000000", 16);
    xLADDER(Q, P, m, A24p, C24);
    point_normalize_coords(Q);

    // x[2^80]P: 98*i + 43
    fp2_print(Q->X, "x[2^80]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 43) && !mpz_cmp_ui(Q->X->b, 98));

    // 2. m = 2^80 - 1
    mpz_init_set_str(m, "ffffffffffffffffffff", 16);
    xLADDER(Q, P, m, A24p, C24);
    point_normalize_coords(Q);

    // x[2^80-1]P: 56*i + 96
    fp2_print(Q->X, "x[2^80-1]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 96) && !mpz_cmp_ui(Q->X->b, 56));

    // 3. m = random in [2^79, 2^80]
    mpz_init_set_str(m, "f5697b000f01c17d4c5e", 16);
    xLADDER(Q, P, m, A24p, C24);
    point_normalize_coords(Q);

    // x[0xf5697b000f01c17d4c5e]P: 94*i + 31
    fp2_print(Q->X, "x[0xf5697b000f01c17d4c5e]P");
    CHECK(!mpz_cmp_ui(Q->X->a, 31) && !mpz_cmp_ui(Q->X->b, 94));

    mpz_clear(m);
}

void test_j_invariant() {
    fp2_t a, c, j_inv;
    fp2_init(&a);
    fp2_init(&c);
    fp2_init(&j_inv);
    fp2_set_uint(c, 1);

    fp2_set_str(a, "92*i + 25");
    fp2_print(a, "a(E)");

    j_invariant(j_inv, a, c);

    // j(E): 9*i + 29
    fp2_print(j_inv, "j(E)");
    CHECK(!mpz_cmp_ui(j_inv->a, 29) && !mpz_cmp_ui(j_inv->b, 9));

    fp2_set_str(a, "125*i + 99");
    fp2_print(a, "a(E)");

    j_invariant(j_inv, a, c);

    // j(E): 79*i + 30
    fp2_print(j_inv, "j(E)");
    CHECK(!mpz_cmp_ui(j_inv->a, 30) && !mpz_cmp_ui(j_inv->b, 79));

    fp2_set_str(a, "43*i + 61");
    fp2_print(a, "a(E)");

    j_invariant(j_inv, a, c);

    // j(E): 78*i + 97
    fp2_print(j_inv, "j(E)");
    CHECK(!mpz_cmp_ui(j_inv->a, 97) && !mpz_cmp_ui(j_inv->b, 78));

    fp2_clear(&a);
    fp2_clear(&c);
    fp2_clear(&j_inv);
}

int main() {
    init_test_variables();

    // Unit tests not related to isogeny - diff with SageMath is not required,
    // therefore tests are run in "SILENT" mode (funciton name is not printed)
    TEST_RUN_SILENT(test_point_set_is_immutable());
    TEST_RUN_SILENT(test_A_A24p_conversion_in_place());
    TEST_RUN_SILENT(test_A_A24p_conversion());

    // p = 431 tests
    set_params_testp431();

    TEST_RUN(test_point_normalize_coords());
    TEST_RUN(test_xDBL_small());
    TEST_RUN(test_xDBLe());
    TEST_RUN(test_xADD_small());
    TEST_RUN(test_xLADDER3PT());

    // p = 139 tests
    set_params_testp139();

    TEST_RUN(test_xLADDER_int());
    TEST_RUN(test_xLADDER());
    TEST_RUN(test_j_invariant());

    clear_test_variables();

    TEST_RUNS_END;
}