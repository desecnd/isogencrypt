#include <stdio.h>

#include "ec_mont.h"
#include "testing.h"

fp2_t A, C;
fp2_t A24_plus, C24;
point_t P, Q;

void init_test_variables() {
    fp2_init(&A);
    fp2_init(&C);
    fp2_init(&A24_plus);
    fp2_init(&C24);
    point_init(&P);
    point_init(&Q);
}

void clear_test_variables() {
    fp2_clear(&A);
    fp2_clear(&C);
    fp2_clear(&A24_plus);
    fp2_clear(&C24);
    point_clear(&P);
    point_clear(&Q);
}

void test_xADD_small() {
    CHECK(!global_fpchar_setup_uint(431));

    // Initialize the Elliptic Curve E: y^2 = x^3 + 6x^2 + x
    fp2_set_uint(A, 6);
    fp2_set_uint(C, 1);
    calc_curve_proj_coeffs(A24_plus, C24, A, C);

    // x(P) = XP/1 = 259 + 271i
    point_set_str(P, "259", "271");
    fp2_print_uint(P->X, "xP");

    // x(Q) = XQ/1 = 262 + 335i
    point_set_str(Q, "262", "335");
    fp2_print_uint(Q->X, "xQ");

    point_t PQdiff, PQsum;
    point_init(&PQdiff);
    point_init(&PQsum);

    // x(P - Q) = xPQdiff = 143 + 411i
    point_set_str(PQdiff, "143", "411");
    fp2_print_uint(PQdiff->X, "xP-Q");

    // Perform the addition:
    xADD(PQsum, P, Q, PQdiff);

    fp2_t x;
    fp2_init(&x);

    // x = X(P+Q)/Z(P+Q) = 61 + 184 * i = x(PQ+)
    fp2_div_unsafe(x, PQsum->X, PQsum->Z);
    fp2_print_uint(x, "xP+Q");

    CHECK(!mpz_cmp_ui(x->a, 416) && !mpz_cmp_ui(x->b, 106));

    fp2_clear(&x);

    point_clear(&PQdiff);
    point_clear(&PQsum);

    CHECK(!global_fpchar_clear());
}

void test_small_xDBL() {

    // Set field characteristic
    CHECK(!global_fpchar_setup_uint(431));

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
    calc_curve_proj_coeffs(A24_plus, C24, A, C);

    fp2_print_uint(A24_plus, "A24p");
    fp2_print_uint(C24, "C24");

    // PX = 292 + 15i, zP = 1
    point_set_str(P, "292", "15");

    fp2_print_uint(P->X, "xP");
    fp2_print_uint(P->Z, "zP");

    // Q = [2]P
    { 
        xDBL(Q, P, A24_plus, C24);
        
        fp2_print_uint(Q->X, "xQ");
        fp2_print_uint(Q->Z, "zQ");

        // Answer in affine coordinates: x(Q) = 61 + 184 * i
        // Which is equal to projective coordinates (this implementation)
        // X(Q) = 157 + 180 * i
        // Z(Q) = 65 + 358 * i

        CHECK(!mpz_cmp_ui(Q->X->a, 157) && !mpz_cmp_ui(Q->X->b, 180));
        CHECK(!mpz_cmp_ui(Q->Z->a, 65) && !mpz_cmp_ui(Q->Z->b, 358));

        fp2_t xQ;
        fp2_init(&xQ);

        // X(Q)/Z(Q) = 61 + 184 * i = x(Q)
        fp2_div_unsafe(xQ, Q->X, Q->Z);
        CHECK(!mpz_cmp_ui(xQ->a, 61) && !mpz_cmp_ui(xQ->b, 184));

        fp2_clear(&xQ);
    }

    CHECK(!global_fpchar_clear());
}

void test_criss_cross_small() {
    CHECK(!global_fpchar_setup_uint(431));
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
    CHECK(!global_fpchar_clear());
}

void test_criss_cross_argsafe() {
    CHECK(!global_fpchar_setup_uint(431));
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
    CHECK(!global_fpchar_clear());
}


int main() {
    init_test_variables();

    TEST_RUN(test_small_xDBL());
    // xADD function
    TEST_RUN(test_xADD_small());

    TEST_RUN(test_criss_cross_small());
    TEST_RUN(test_criss_cross_argsafe());
    TEST_RUNS_END;

    clear_test_variables();
    return 0;
}